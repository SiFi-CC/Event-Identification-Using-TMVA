#include "EIComp.hh"
#include <iostream>
using namespace std;

int main(int argc, char* argv[]){
  //Int_t nev = 100000;
  //Int_t nev = 2200797;
   //Int_t nev = 1100398; // half number
  //Int_t nev = 1278098;
  //Int_t nev = 639049; // half number
  //Int_t nev = 526474;  
  
  if (argc != 2) {
    cout << "To run type: ./CompId path_to_EI.txt" << endl;
    return 0;
  }
  
  TString path(argv[1]);
  
  EIComp* id;
  
  try {
      
    id = new EIComp(path);
  } catch (const char* message) {
      cout << message << endl;
      return 0;
  }
 
  id->Identify();
  
  
  
  delete id;
 
  return 1;
}
