#include "TrueDaughter.h"

TrueDaughter::TrueDaughter(){

  j_IsEmpty = true;
  j_PDG = -9999;
  j_ID = -9999;
  j_startX = -9999.;
  j_startY = -9999.;
  j_startZ = -9999.;
  j_startPx = -9999.;
  j_startPy = -9999.;
  j_startPz = -9999.;
  j_startP = -9999.;
  j_endX = -9999.;
  j_endY = -9999.;
  j_endZ = -9999.;
  j_Process = "";
  j_endProcess = "";
}

TrueDaughter::~TrueDaughter(){}

void TrueDaughter::Set_IsEmpty(bool i_IsEmpty){ j_IsEmpty = i_IsEmpty; }
void TrueDaughter::Set_PDG(int i_PDG){ j_PDG = i_PDG; }
void TrueDaughter::Set_ID(int i_ID){ j_ID = i_ID; }
void TrueDaughter::Set_startX(double i_startX){ j_startX = i_startX; }
void TrueDaughter::Set_startY(double i_startY){ j_startY = i_startY; }
void TrueDaughter::Set_startZ(double i_startZ){ j_startZ = i_startZ; }
void TrueDaughter::Set_startPx(double i_startPx){ j_startPx = i_startPx; }
void TrueDaughter::Set_startPy(double i_startPy){ j_startPy = i_startPy; }
void TrueDaughter::Set_startPz(double i_startPz){ j_startPz = i_startPz; }
void TrueDaughter::Set_startP(double i_startP){ j_startP = i_startP; }
void TrueDaughter::Set_endX(double i_endX){ j_endX = i_endX; }
void TrueDaughter::Set_endY(double i_endY){ j_endY = i_endY; }
void TrueDaughter::Set_endZ(double i_endZ){ j_endZ = i_endZ; }
void TrueDaughter::Set_Process(string i_Process){ j_Process = i_Process; }
void TrueDaughter::Set_endProcess(string i_endProcess){ j_endProcess = i_endProcess; }
