#ifndef TrueDaughter_h
#define TrueDaughter_h

#include <string>
#include <vector>

using namespace std;

class TrueDaughter{
public:

  TrueDaughter();
  virtual ~TrueDaughter();

  // == Empty
  void Set_IsEmpty(bool i_IsEmpty);
  inline bool IsEmpty() const { return j_IsEmpty; }

  // == Variables
  void Set_PDG(int i_PDG);
  void Set_ID(int i_ID);
  void Set_len(int i_len);
  void Set_startX(double i_startX);
  void Set_startY(double i_startY);
  void Set_startZ(double i_startZ);
  void Set_startPx(double i_startPx);
  void Set_startPy(double i_startPy);
  void Set_startPz(double i_startPz);
  void Set_startP(double i_startP);
  void Set_endX(double i_endX);
  void Set_endY(double i_endY);
  void Set_endZ(double i_endZ);
  void Set_Process(string i_Process);
  void Set_endProcess(string i_endProcess);

  inline int PDG() const { return j_PDG; }
  inline int ID() const { return j_ID; }
  inline double startX() const { return j_startX; }
  inline double	startY() const { return j_startY; }
  inline double	startZ() const { return j_startZ; }
  inline double	startPx() const { return j_startPx; }
  inline double startPy() const { return j_startPy; }
  inline double startPz() const { return j_startPz; }
  inline double startP() const { return j_startP; }
  inline double	endX() const { return j_endX; }
  inline double endY() const { return j_endY; }
  inline double endZ() const { return j_endZ; }
  inline string Process() const { return j_Process; }
  inline string endProcess() const { return j_endProcess; }
  
private:
  bool j_IsEmpty;
  int j_PDG;
  int j_ID;
  double j_startX;
  double j_startY;
  double j_startZ;
  double j_startPx;
  double j_startPy;
  double j_startPz;
  double j_startP;
  double j_endX;
  double j_endY;
  double j_endZ;
  string j_Process;
  string j_endProcess;  
};

#endif
