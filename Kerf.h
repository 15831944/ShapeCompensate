#ifndef SHAPECOMPENSATE_KERF_H__
#define SHAPECOMPENSATE_KERF_H__

#include <vector>

#include "GCodeStruct.h"

class Kerf {
 public:
  Kerf();
  ~Kerf();
  
  void SetKerfValue(double kerf_value);

  bool ReadGCode(const std::string &file_name, std::vector<std::string> &code_lines);
  bool WriteGCode(const std::string &file_name, const std::vector<std::string> &contents);

  void GCodeParse(const std::vector<std::string> &code_lines);
  void GenerateKerfGCode(std::vector<std::string> &code_lines);

 private:
  bool is_absolute;
  double kerf;
  std::vector<GCodeARRAY_STRUCT> GfileFloatKerf; // 有割缝的浮点型切割代码
  std::vector<GCodeARRAY_STRUCT> GfileFloatNoKerf; // 没有割缝的浮点型切割代码
  GraphyLimit graphylimitxy;

  void CircleCheFen(std::vector<GCodeARRAY_STRUCT> &GCodeArry);
  int IgnoreLittleLine(std::vector<GCodeARRAY_STRUCT> &GCodeArry);
  double GetRadius(GCodeARRAY_STRUCT &pG);
  double GetTangent(GCodeARRAY_STRUCT *pG, int StartOrEnd);
  int GetAddKerfGCode(GCodeARRAY_STRUCT &pNoKerfG,
      GCodeARRAY_STRUCT *AddKerfGCode, double kerfvalue, int dir);

  int AddOrTrunc(GCodeARRAY_STRUCT &pPreviousLine, GCodeARRAY_STRUCT &pNextLine,
      GCodeARRAY_STRUCT *pAddLine, double kerfvalue, int dir);

  int Setupkerf(GCodeARRAY_STRUCT *pGcode, double *dx, double *dy,
      double kerfvlaue, int dir);

  int Canclekerf(GCodeARRAY_STRUCT *pGcode, double *dx, double *dy,
      double kerfvlaue, int dir);

  double GetCodeValue(const std::string &code_line,
                      const std::string &match);

  std::string TransformGCodeLine(const GCodeARRAY_STRUCT &gcode_array);

  void g2kerf(std::vector<GCodeARRAY_STRUCT> &DesKerFile,
              std::vector<GCodeARRAY_STRUCT> &NoKerfFile);

}; // class Kerf

#endif // SHAPECOMPENSATE_KERF_H__
