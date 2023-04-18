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
struct Coordinate {
    double x, y, z;
};
void replace_coordinates_in_pdb(const string &pdbFileName, const string &coordFileName, const string
&outputFileName) {
    ifstream pdbFile(pdbFileName);
    ifstream coordFile(coordFileName);
    ofstream outputFile(outputFileName);

    if (!pdbFile.is_open() || !coordFile.is_open()) {
        cerr << "Error opening input files." << endl;
        return;
    }

    vector<Coordinate> newCoords;
    double x, y, z;
    while (coordFile >> x >> y >> z) {
        newCoords.push_back({x, y, z});
    }

    string line;
    size_t index = 0;

    while (getline(pdbFile, line)) {
        if (index < newCoords.size() && line.substr(0, 4) == "ATOM") {
            string beginning = line.substr(0, 30);
            string ending = line.substr(54);
            outputFile << beginning
                       << fixed << setprecision(3) << setw(8) << newCoords[index].x
                       << setw(8) << newCoords[index].y << setw(8) << newCoords[index].z
                       << ending << endl;
            index++;
        } else {
            outputFile << line << endl;
        }
    }

    pdbFile.close();
    coordFile.close();
    outputFile.close();
}




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

    wholeModel.readPDBfile(fileModel, PDB::AllSelector());
    wholeTarget.readPDBfile(fileTarget, PDB::AllSelector());
    fileModel.clear();
    fileTarget.clear();
    fileModel.seekg(0, fileModel.beg);
    fileTarget.seekg(0, fileTarget.beg);
    modelBackBone.readPDBfile(fileModel, PDB::CAlphaSelector());
    targetBackBone.readPDBfile(fileTarget, PDB::CAlphaSelector());
    //TODO: choose phosphate for rna here

    // calculate center of mass
    //TODO: DO THIS FOR JUST ONE AND DROP THE SECOND ON IT
    Vector3 vectTargetMass(0, 0, 0);
    for (unsigned int i = 0; i < targetBackBone.size(); i++) {
        vectTargetMass += targetBackBone[i].position();
    }
    vectTargetMass /= targetBackBone.size();

    // transform the molecules to the center of the coordinate system
    targetBackBone += (-vectTargetMass);

    // next we insert the target backcone into hash
    // this will help us to find atoms that are close faster
    GeomHash<Vector3, int> gHash(3, m_fDistThr); // 3 is a dimension and m_fDistThr is the size of the hash cube
    for (unsigned int i = 0; i < targetBackBone.size(); i++) {
        gHash.insert(targetBackBone[i].position(), i); // coordinate is the key to the hash, we store atom index
    }
    cout << "Size of targetBackBone: " << targetBackBone.size() << endl;
//	for (unsigned int i = 0; i < wholeTarget.size(); i++) {
//		cout << "Atom " << i << " position in wholeTarget: " << wholeTarget[i].position() << endl;
//		// Try to retrieve the atom position from gHash using the query method
//		HashResult<int> result;
//		gHash.query(wholeTarget[i].position(), m_fDistThr, result);
//
//		if (!result.empty()) {
//			int atomIndex = *result.begin();
//			cout << "Atom " << atomIndex << " position in gHash: " << wholeTarget[atomIndex].position() << endl;
//		} else {
//			cout << "Atom " << i << " not found in gHash" << endl;
//		}
//	}


    unsigned int iMaxSize = 0;
    RigidTrans3 rtransBest;
    Match bestMatch;
    Match match;

    for (unsigned int i = 0; i < targetBackBone.size() - 2; i++) {
        Vector3 a_t = targetBackBone[i].position();
        Vector3 b_t = targetBackBone[i + 1].position();
        Vector3 c_t = targetBackBone[i + 2].position();

        for (unsigned int j = 0; j < modelBackBone.size() - 2; j++) {
            Vector3 a_m = modelBackBone[j].position();
            Vector3 b_m = modelBackBone[j + 1].position();
            Vector3 c_m = modelBackBone[j + 2].position();

            RigidTrans3 trans_target = RigidTrans3(a_t, b_t, c_t);
            RigidTrans3 trans_model = RigidTrans3(a_m, b_m, c_m);
            const RigidTrans3 final_trans = (!trans_target) * trans_model;

            // apply the transformation to all model backbone molecules
            for (unsigned int k = 0; k < modelBackBone.size(); k++) {
                Vector3 mol = modelBackBone[k].position();
                mol = final_trans * mol;

                // find close target molecule atoms using the hash
                HashResult<int> result;
                gHash.query(mol, m_fDistThr, result); // key is mol atom coordinate

                // Print out the result
//				cout << "Close target molecule atoms for model atom " << k << ":" << endl;
//				for (auto x = result.begin(); x != result.end(); x++) {
//					cout << "is  Target atom index: " << *x << endl;
//				}

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
//                result.clear();
            }

            // the match class will now check for any pair clash and optimize the fit
            match.calculateBestFit(targetBackBone, modelBackBone);

//			cout << "Rigid Trans: " << match.rigidTrans() << '\n';

//			int matchingParticle;
//			for(int m=0; m<modelBackBone.size(); m++){
//				if ((matchingParticle = match.sceneParticle(m)) >= 0){
//					cout << m << '\t' << matchingParticle << '\n';
//				}
//			}

            if (iMaxSize < match.size()) {
                iMaxSize = match.size();
                rtransBest = match.rigidTrans();
                bestMatch = match;
            }
        }
    }

    cout << "Max Alignment Size: " << iMaxSize << endl;
    cout << "model size: " << modelBackBone.size() << endl;
    cout << "target size: " << targetBackBone.size() << endl;
    cout << bestMatch.size() << '\t' << bestMatch.calculateTotalScore() << '\t' << rtransBest << '\t' << endl;

    auto end = chrono::system_clock::now();

    chrono::duration<double> elapsed_seconds = end - start;
    cout << "elapsed time: " << elapsed_seconds.count() << "s" << endl;

    fileModel.close();
    fileTarget.close();

    // apply the transformation (now not in a loop) to all model and target molecules and export them to pdb's
    wholeModel *= rtransBest;
    wholeTarget += (-vectTargetMass);

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


    ////////////////////////////////////// PDB WRITING //////////////////////////////////////
    replace_coordinates_in_pdb("input_files/1m5oC.pdb", "model.txt", "1m5oC-model_modified.pdb");
    replace_coordinates_in_pdb("input_files/1b7fA.pdb", "target.txt", "1b7fA-target_modified.pdb");

}