// Copyright 2016 Fangling Software Co., Ltd. All Rights Reserved.
// Author: shizhan-shelly@hotmail.com (Zhan Shi)

#ifndef SHAPECOMPENSATE_KERF_H__
#define SHAPECOMPENSATE_KERF_H__

#include <vector>

#include "GCodeDefinition.h"

class Kerf {
 public:
  Kerf();
  ~Kerf();

  void SetKerfValue(double kerf_value);

  void GKerfProc(const std::string &noKerfFile, const std::string &KerfFile);

 private:
  double param_kerf;
  std::vector<GCodeStruct> GfileFloatKerf; // �и��ĸ������и����
  std::vector<GCodeStruct> GfileFloatNoKerf; // û�и��ĸ������и����
  GraphyLimit graphylimitxy;

  void CircleCheFen(std::vector<GCodeStruct> &GCodeArry);
  int IgnoreLittleLine(std::vector<GCodeStruct> &GCodeArry);
  double GetRadius(GCodeStruct &pG);
  double GetTangent(GCodeStruct &pG, int StartOrEnd);
  int GetAddKerfGCode(GCodeStruct &pNoKerfG,
      GCodeStruct &AddKerfGCode, double kerfvalue, int dir);

  int AddOrTrunc(GCodeStruct &pPreviousLine, GCodeStruct &pNextLine,
      GCodeStruct &pAddLine, int dir);

  int Setupkerf(GCodeStruct &pGcode, double &dx, double &dy,
      double kerfvlaue, int dir);

  int Canclekerf(GCodeStruct &pGcode, double &dx, double &dy,
      double kerfvlaue, int dir);

  void g2kerf(std::vector<GCodeStruct> &DesKerFile,
              std::vector<GCodeStruct> &NoKerfFile);

}; // class Kerf

#endif // SHAPECOMPENSATE_KERF_H__
