#pragma once

#include <vector>
#include <memory>


class cProtein;
class cWaterResidue;

std::vector<std::unique_ptr<cWaterResidue>> GetHydrationShell(cProtein *protein,
                                                              double sa_radius,
                                                              double sa_smoothed_margin,
                                                              double shell_width);
