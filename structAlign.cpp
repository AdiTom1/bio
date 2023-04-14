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

int main(int argc , char* argv[]){
	// measure the run time
	auto start = std::chrono::system_clock::now();

	if(argc !=4) {
		std::cerr << "Usage: "<<argv[0]<< "structalign epsilon pdb1 pdb2" << std::endl;
		exit(1);
	}

	//********Parameters********************
	float m_fDistThr = atoi(argv[1]); // number of random rotations to try

	std::cout << "Distance threshold: "<< m_fDistThr  << std::endl;

	// read the two files into Molecule
	Molecule<Atom> molModel, molTarget;

	std::ifstream fileModel(argv[3]);
	std::ifstream fileTarget(argv[2]);

	if(!fileModel) {
		std::cout<< "File " << argv[3] << "does not exist." << std::endl;
		return 0;
	}
	if(!fileTarget) {
		std::cout << "File " << argv[2] << "does not exist." << std::endl;
		return 0;
	}

	molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
	molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());
	//TODO: choose phosphate for rna here

	// calculate center of mass
	//TODO: DO THIS FOR JUST ONE AND DROP THE SECOND ON IT
//	Vector3 vectModelMass(0,0,0);
//	for(unsigned int i=0; i<molModel.size(); i++) {
//		vectModelMass+=molModel[i].position();
//	}
//	vectModelMass/=molModel.size();
	Vector3 vectModelMass(0,0,0);
	for(unsigned int i=0; i<molModel.size(); i++) {
		vectModelMass+=molModel[i].position();
	}
	vectModelMass/=molModel.size();


	Vector3 vectTargetMass(0,0,0);
	for(unsigned int i=0; i<molTarget.size(); i++) {
		vectTargetMass+=molTarget[i].position();
	}
	vectTargetMass/=molTarget.size();

	// transform the molecules to the center of the coordinate system
//	molModel+=(-vectModelMass);
	molTarget+=(-vectTargetMass);


	// next we insert the target molecule into hash
	// this will help us to find atoms that are close faster
	GeomHash <Vector3,int> gHash(3,m_fDistThr); // 3 is a dimension and m_fDistThr is the size of the hash cube
	for(unsigned int i=0; i<molTarget.size(); i++) {
		gHash.insert(molTarget[i].position(), i); // coordinate is the key to the hash, we store atom index
	}


	// now we try random rotations and choose the best alignment from random rotations
	unsigned int iMaxSize=0;
	RigidTrans3 rtransBest;

	for(unsigned int i=0; i< molTarget.size()-2; i++) {
		Vector3 a_t = molTarget[i].position();
		Vector3 b_t = molTarget[i+1].position();
		Vector3 c_t = molTarget[i+2].position();
//		Triangle targetTriangle = Triangle(a_t,b_t,c_t);
		// match is a class that stores the correspondence list, eg.
		// pairs of atoms, one from each molecule, that are matching

		// apply rotation on each atom in the model molecule and
		// add the pairs of atoms (one from target and one from model)
		// that are close enough to the match list
		for(unsigned int j=0; j< molModel.size()-2; j++)
		{

			Vector3 a_m = molTarget[j].position();
			Vector3 b_m = molTarget[j + 1].position();
			Vector3 c_m = molTarget[j + 2].position();
//			Triangle modelTriangle = Triangle(a_m,b_m,c_m);

			RigidTrans3 trans_target = RigidTrans3(a_t, b_t, c_t);
			RigidTrans3 trans_model = RigidTrans3(a_m, b_m, c_m);
			const RigidTrans3 final_trans = (!trans_target) * trans_model;


			Match match;
			for (unsigned int k = 0; k < molModel.size(); k++)
			{
				Vector3 mol = molModel[k].position();
				mol *= final_trans;

				// find close target molecule atoms using the hash
				HashResult<int> result;
				gHash.query(mol, m_fDistThr, result); // key is mol atom coordinate

				// check if the atoms in the result are inside the distance threshold
				// the hash is a cube shape, there can be atoms further that the threshold
				for (auto x = result.begin(); x != result.end(); x++)
				{
					float dist = mol.dist(molTarget[*x].position());
					if (dist <= m_fDistThr)
					{
						float score = (1 / (1 + dist));
						match.add(*x, i, score, score);
					}
				}
				result.clear();
			}


			//calculates transformation that is a little better than "rotation"
			match.calculateBestFit(molTarget, molModel);

			if (iMaxSize < match.size())
			{
				iMaxSize = match.size();
				rtransBest = match.rigidTrans();
			}
		}
	}

	std::cout << "Max Alignment Size: " << iMaxSize << std::endl;
	std::cout << "Rigid Trans: " <<
			  RigidTrans3(Vector3(0,0,0),vectTargetMass)*
			  rtransBest*
			  RigidTrans3(Vector3(0,0,0),(-vectModelMass)) << std::endl;

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;

}
