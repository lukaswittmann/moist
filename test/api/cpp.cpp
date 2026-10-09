#include "moist.h"
#include <cstdio>
#include <type_traits>

static_assert(std::is_same<decltype(moist_check_error(nullptr)), int>::value,
              "moist_check_error must return int");
static_assert(moist_success == 0 && moist_failure == 1 && moist_invalid_error == 2,
              "status values are ABI");

template<class Handle>
bool delete_handle(Handle& handle, const char* name)
{
    if (!handle) {
        std::fprintf(stderr, "%s constructor returned null\n", name);
        return false;
    }
    moist_delete(handle);
    if (handle) {
        std::fprintf(stderr, "%s delete did not clear handle\n", name);
        return false;
    }
    moist_delete(handle);
    return handle == nullptr;
}

int fail(int line)
{
    std::fprintf(stderr, "C++ API check failed at line %d\n", line);
    return 1;
}

int main()
{
    moist_error error = moist_new_error();
    if (!error || moist_check_error(error) != moist_success ||
        moist_check_error(nullptr) != moist_invalid_error) return fail(__LINE__);
    moist_context context = moist_new_context(error, 2, 0, false);
    if (!context || moist_check_error(error) != moist_success) return fail(__LINE__);
    moist_lsf lsf = moist_new_svdw_lsf(error, nullptr);
    if (!lsf || moist_check_error(error) != moist_success) return fail(__LINE__);
    moist_cavity cavity = moist_new_drop_cavity(error, context, lsf, nullptr, nullptr);
    if (!cavity || moist_check_error(error) != moist_success) return fail(__LINE__);
    const int numbers[] = {1};
    const double positions[] = {0.0, 0.0, 0.0};
    moist_structure structure = moist_new_structure(error, 1, numbers, positions,
                                                   nullptr, nullptr);
    if (!structure || moist_check_error(error) != moist_success) return fail(__LINE__);
    moist_radii radii = moist_new_cpcm_radii(error);
    if (!radii || moist_check_error(error) != moist_success) return fail(__LINE__);
    moist_model model = moist_new_model(error, context, cavity);
    if (!model || moist_check_error(error) != moist_success) return fail(__LINE__);
    // Cavity and model retain the context past its handle
    if (!delete_handle(context, "context")) return fail(__LINE__);
    moist_component component = moist_new_pv_component(error, nullptr, 0.001);
    if (!component || moist_check_error(error) != moist_success) return fail(__LINE__);
    moist_add_model_component(error, model, component);
    if (moist_check_error(error) != moist_success) return fail(__LINE__);
    moist_update_model(error, model, structure);
    if (moist_check_error(error) != moist_success) return fail(__LINE__);
    moist_coupling coupling = moist_new_coupling(error, model);
    if (!coupling || moist_check_error(error) != moist_success) return fail(__LINE__);
    moist_response response = moist_new_response(error);
    if (!response || moist_check_error(error) != moist_success) return fail(__LINE__);
    if (!delete_handle(coupling, "coupling") ||
        !delete_handle(response, "response") ||
        !delete_handle(component, "component") ||
        !delete_handle(model, "model") ||
        !delete_handle(structure, "structure") ||
        !delete_handle(radii, "radii") ||
        !delete_handle(lsf, "lsf") ||
        !delete_handle(cavity, "cavity")) return fail(__LINE__);
    moist_cavity failed = moist_new_drop_cavity(error, nullptr, nullptr, nullptr, nullptr);
    if (moist_check_error(error) != moist_failure || failed != nullptr) return fail(__LINE__);
    return delete_handle(error, "error") ? 0 : 1;
}
