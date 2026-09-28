#include "moist.h"
#include <type_traits>

static_assert(std::is_same<decltype(moist_check_error(nullptr)), int>::value,
              "moist_check_error must return int");

int main()
{
    moist_error error = moist_new_error();
    moist_lsf lsf = moist_new_svdw_lsf(error, nullptr);
    moist_cavity cavity = moist_new_drop_cavity(error, lsf, nullptr, nullptr);
    const int status = moist_check_error(error);
    moist_delete(lsf);
    moist_delete(cavity);
    moist_structure structure = nullptr;
    moist_radii radii = nullptr;
    moist_model model = nullptr;
    moist_component component = nullptr;
    moist_coupling coupling = nullptr;
    moist_response response = nullptr;
    moist_delete(structure);
    moist_delete(radii);
    moist_delete(model);
    moist_delete(component);
    moist_delete(coupling);
    moist_delete(response);
    moist_delete(error);
    moist_delete(error);
    return status || error || cavity || lsf;
}
