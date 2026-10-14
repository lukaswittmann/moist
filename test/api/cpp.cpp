#include "moist.h"
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

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

// Host timer stopped when its scope ends, on every path out of it (C++ Core
// Guidelines R.1, E.19). A refused start throws std::runtime_error with
// MOIST's message. The stop reports into its own moist_error, so the caller's
// diagnostic survives, and the destructor never throws (C.36, C.37). The
// context must outlive the guard.
class host_timer {
public:
    host_timer(moist_context context, std::string name)
        : context_(context), name_(std::move(name))
    {
        moist_error error = moist_new_error();
        moist_start_context_timer(error, context_, name_.c_str());
        const bool started = moist_check_error(error) == moist_success;
        char message[512] = {0};
        const int capacity = static_cast<int>(sizeof message);
        if (!started) moist_get_error(error, message, &capacity);
        moist_delete(error);
        if (!started) throw std::runtime_error(message);
    }
    ~host_timer()
    {
        moist_error diagnostic = moist_new_error();
        moist_stop_context_timer(diagnostic, context_, name_.c_str());
        if (moist_check_error(diagnostic) != moist_success) {
            char message[512] = {0};
            const int capacity = static_cast<int>(sizeof message);
            moist_get_error(diagnostic, message, &capacity);
            std::fprintf(stderr, "%s\n", message);
        }
        moist_delete(diagnostic);
    }
    host_timer(const host_timer&) = delete;
    host_timer& operator=(const host_timer&) = delete;

private:
    moist_context context_;
    std::string name_;
};

// Whether the diagnostic in `error` contains `text`
bool error_mentions(moist_error error, const char* text)
{
    char message[512] = {0};
    const int capacity = static_cast<int>(sizeof message);
    if (moist_check_error(error) == moist_success) return false;
    moist_get_error(error, message, &capacity);
    return std::strstr(message, text) != nullptr;
}

// The guard stops its timer when the scope throws, keeps the caller's
// diagnostic, and refuses a name breaking the naming rules
bool host_timer_stops_on_every_path(moist_error error, moist_context context)
{
    int count = 0;
    try {
        host_timer timer(context, "Fock build");
        moist_stop_context_timer(error, context, "other");
        if (!error_mentions(error, "'Fock build'")) return false;
        moist_get_context_num_threads(error, nullptr, &count);
        throw std::logic_error("host work failed");
    } catch (const std::logic_error&) {
    }
    if (!error_mentions(error, "Context handle is missing")) return false;
    moist_stop_context_timer(error, context, "Fock build");
    if (!error_mentions(error, "No host timer is open")) return false;
    try {
        host_timer timer(context, "A ");
        return false;
    } catch (const std::runtime_error& refused) {
        return std::strstr(refused.what(), "space") != nullptr;
    }
}

int main()
{
    moist_error error = moist_new_error();
    if (!error || moist_check_error(error) != moist_success ||
        moist_check_error(nullptr) != moist_invalid_error) return fail(__LINE__);
    moist_context context = moist_new_context(error, 2, 0, false);
    if (!context || moist_check_error(error) != moist_success) return fail(__LINE__);
    if (!host_timer_stops_on_every_path(error, context)) return fail(__LINE__);
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
    // The build phase of a PV model on a geometric cavity asks for nothing
    moist_prepare_model_build(error, model, coupling);
    if (moist_check_error(error) != moist_success) return fail(__LINE__);
    if (moist_next_coupling_request(error, coupling)) return fail(__LINE__);
    moist_build_model(error, model, coupling);
    if (moist_check_error(error) != moist_success) return fail(__LINE__);
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
