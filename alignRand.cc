#include "Vector3.h"
#include "Atom.h"
#include "RigidTrans3.h"
#include "Matrix3.h"
#include "Molecule.h"
#include "PDB.h"
#include "Match.h"
#include "GeomHash.h"
#include <chrono>
#include <iostream>

int main(int argc , char* argv[]){
  // measure the run time
  auto start = std::chrono::system_clock::now();
  
  if(argc !=5) {
    std::cerr << "Usage: "<<argv[0]<< " target.pdb model.pdb num_rotations dist_threshold" << std::endl;
    exit(1);
  }

  //********Parameters********************
  int m_iRotations = atoi(argv[3]); // number of random rotations to try
  float m_fDistThr = atof(argv[4]); // distance threshold on atoms in correspondence

  std::cout << "Number of rotations: "<< m_iRotations << std::endl;
  std::cout << "Distance threshold: "<< m_fDistThr  << std::endl;

  // read the two files into Molecule
  Molecule<Atom> molModel, molTarget;

  std::ifstream fileModel(argv[2]);
  std::ifstream fileTarget(argv[1]);

  if(!fileModel) {
    std::cout<< "File " << argv[1] << "does not exist." << std::endl;
    return 0;
  }
  if(!fileTarget) {
    std::cout << "File " << argv[2] << "does not exist." << std::endl;
    return 0;
  }

  molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
  molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());

  // calculate center of mass
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
  molModel+=(-vectModelMass);
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
  for(int iRand=0; iRand < m_iRotations; iRand++) {
    // random rotation matrix
    Matrix3 rotation((drand48()-0.5)*2*3.1415,
                     (drand48()-0.5)*2*3.1415,
                     (drand48()-0.5)*2*3.1415); 
    
    // match is a class that stores the correspondence list, eg.
    // pairs of atoms, one from each molecule, that are matching 
    Match match;

    // apply rotation on each atom in the model molecule and
    // add the pairs of atoms (one from target and one from model) 
    // that are close enough to the match list
    for(unsigned int i=0; i< molModel.size(); i++) {
      Vector3 mol_atom = rotation*molModel[i].position(); // rotate

      // find close target molecule atoms using the hash
      HashResult<int> result; 
      gHash.query(mol_atom, m_fDistThr, result); // key is mol atom coordinate

      // check if the atoms in the result are inside the distance threshold
      // the hash is a cube shape, there can be atoms further that the threshold 
      for(auto x = result.begin(); x != result.end(); x++) {
        float dist = mol_atom.dist(molTarget[*x].position());
	    if(dist <= m_fDistThr) {
          float score = (1 / (1 + dist));
          match.add( *x , i, score, score );
        }
      }
      result.clear();
    }

    //calculates transformation that is a little better than "rotation"
    match.calculateBestFit(molTarget, molModel);

    if(iMaxSize < match.size() ){
      iMaxSize = match.size();
      rtransBest=match.rigidTrans();
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
