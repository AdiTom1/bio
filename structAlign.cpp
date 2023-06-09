//
// Created by aloni on 4/17/2023.
//

#include "Vector3.h"
#include "Atom.h"
#include "RigidTrans3.h"
#include "Matrix3.h"
#include "Molecule.h"
#include "Triangle.h"
#include "PDB.h"
#include "Match.h"
#include "GeomHash.h"
#include <chrono>
#include <iostream>
using namespace std;


// for exporting to pdb
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>



int main(int argc, char* argv[]) {
    // measure the run time
    auto start = chrono::system_clock::now();

    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " structalign epsilon pdb1 pdb2" << endl;
        exit(1);
    }

    //********Parameters********************
    float m_fDistThr = strtof(argv[1], nullptr); // distance threshold

    cout << "Distance threshold: " << m_fDistThr << endl;

    // read the two files into Molecule
    Molecule<Atom> wholeModel, wholeTarget, modelBackBone, targetBackBone;

    ifstream fileTarget(argv[2]);
    ifstream fileModel(argv[3]);

    if (!fileModel) {
        cout << "File " << argv[3] << " does not exist." << endl;
        return 0;
    }
    if (!fileTarget) {
        cout << "File " << argv[2] << " does not exist." << endl;
        return 0;
    }
    string first_line;
    bool isRNA = false;
    fileModel.clear();
    fileTarget.clear();
    fileModel.seekg(0, fileModel.beg);
    fileTarget.seekg(0, fileTarget.beg);

    getline(fileModel, first_line);
    if (first_line.find("RNA") != string::npos){
        isRNA = true;
    }

    wholeModel.readPDBfile(fileModel, PDB::AllSelector());
    wholeTarget.readPDBfile(fileTarget, PDB::AllSelector());
    fileModel.clear();
    fileTarget.clear();
    fileModel.seekg(0, fileModel.beg);
    fileTarget.seekg(0, fileTarget.beg);

    if (isRNA){
        modelBackBone.readPDBfile(fileModel, PDB::PSelector());
        targetBackBone.readPDBfile(fileTarget, PDB::PSelector());
    }else{
            modelBackBone.readPDBfile(fileModel, PDB::CAlphaSelector());
            targetBackBone.readPDBfile(fileTarget, PDB::CAlphaSelector());

    };



    // next we insert the target backcone into hash
    // this will help us to find atoms that are close faster
    GeomHash<Vector3, int> gHash(3, m_fDistThr); // 3 is a dimension and m_fDistThr is the size of the hash cube
    for (unsigned int i = 0; i < targetBackBone.size(); i++) {
        gHash.insert(targetBackBone[i].position(), i); // coordinate is the key to the hash, we store atom index
    }
    cout << "Size of targetBackBone: " << targetBackBone.size() << endl;


    unsigned int iMaxSize = 0;
    RigidTrans3 rtransBest;
    Match bestMatch;


    for (unsigned int i = 0; i < targetBackBone.size() - 2; i++) {

        Vector3 a_t = targetBackBone[i].position();
        Vector3 b_t = targetBackBone[i + 1].position();
        Vector3 c_t = targetBackBone[i + 2].position();
        Triangle targetTriangle= Triangle(a_t, b_t, c_t);

        for (unsigned int j = 0; j < modelBackBone.size() - 2; j++) {
            Match match;
            Vector3 a_m = modelBackBone[j].position();
            Vector3 b_m = modelBackBone[j + 1].position();
            Vector3 c_m = modelBackBone[j + 2].position();
            Triangle modelTriangle= Triangle(a_m, b_m, c_m);

            const RigidTrans3 final_trans =   targetTriangle | modelTriangle ;

            // apply the transformation to all model backbone molecules
            for (unsigned int k = 0; k < modelBackBone.size(); k++) {
                Vector3 mol = modelBackBone[k].position();
                mol = final_trans * mol;

                // find close target molecule atoms using the hash
                HashResult<int> result;
                gHash.query(mol, m_fDistThr, result); // key is mol atom coordinate

                // check if the atoms in the result are really inside the distance threshold
                // the hash is cube shaped, there can be atoms further than the threshold.
                // choose the closest one add
                float dist;
                for (auto x = result.begin(); x != result.end(); x++) {
                    dist = mol.dist(targetBackBone[*x].position());
                    if (dist <= m_fDistThr) {
                        float score = (1 / (1 + dist));
                        match.add(*x, k, score, score);
                    }
                }
                result.clear();
            }

            // the match class will now check for any pair clash and optimize the fit
            match.calculateBestFit(targetBackBone, modelBackBone);

            if (iMaxSize < match.size()) {
                iMaxSize = match.size();
                bestMatch = match;
            }
        }
    }

    cout << "Max Alignment Size: " << iMaxSize << endl;
    cout << "model size: " << modelBackBone.size() << endl;
    cout << "target size: " << targetBackBone.size() << endl;
    cout << bestMatch.size() << '\t' << bestMatch.rmsd() << '\t' << bestMatch.rigidTrans() << '\t' << endl;

    auto end = chrono::system_clock::now();

    chrono::duration<double> elapsed_seconds = end - start;
    cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;

    fileModel.close();
    fileTarget.close();

    // apply the transformation (now not in a loop) to all model and target molecules and export them to pdb's
    wholeModel *= bestMatch.rigidTrans();

    ofstream model;
    model.open ("model.txt");
    for (unsigned int i = 0; i < wholeModel.size(); i++) {
        model << wholeModel[i].position() << endl;
    }
    model.close();

    ofstream target;
    target.open ("target.txt");
    for (unsigned int i = 0; i < wholeTarget.size(); i++) {
        target << wholeTarget[i].position() << endl;
    }
    target.close();


    std::ofstream output_file( "transformed.pdb");
    if(output_file.is_open()){
        output_file << wholeModel;
        output_file.close();
    }

}