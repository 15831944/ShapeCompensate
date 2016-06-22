#include "Kerf.h"

int main(int argc, char *argv[]) {
  Kerf kerf;
  kerf.SetKerfValue(1.5);
  std::vector<std::string> lines;
  kerf.ReadGCode("d://MRT20GL19.TXT", lines);

  kerf.GCodeParse(lines);

  std::vector<std::string> new_lines;
  kerf.GenerateKerfGCode(new_lines);
  kerf.WriteGCode("d://cutfile.ker", new_lines);
  return 1;
}
