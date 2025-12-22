#ifndef FIELDS_COMPOSITE_FIELD_H_
#define FIELDS_COMPOSITE_FIELD_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/pde.h"
#include "core/operation/operation.h"

namespace bdm {
namespace skibidy {

// Container that owns an array of PDE objects.
// Provides named lookup, batch initialization, and per-step source dispatch.
class CompositeField {
 public:
  void Add(std::unique_ptr<PDE> pde) {
    const char* name = pde->GetName();
    index_[name] = pdes_.size();
    pdes_.push_back(std::move(pde));
  }

  PDE* Get(const std::string& name) const {
    auto it = index_.find(name);
    if (it == index_.end()) return nullptr;
    return pdes_[it->second].get();
  }

  DiffusionGrid* Grid(const std::string& name, Simulation* sim) const {
    auto* pde = Get(name);
    return pde ? pde->Grid(sim) : nullptr;
  }

  void InitAll(Simulation* sim) {
    for (auto& p : pdes_) {
      p->Init(sim);
    }
  }

  void ApplyAllSources(Simulation* sim) {
    for (auto& p : pdes_) {
      p->ApplySource(sim, *this);
    }
  }

  void ApplyWoundAll(Simulation* sim, real_t cx, real_t cy, real_t r) {
    for (auto& p : pdes_) {
      p->ApplyWound(sim, cx, cy, r);
    }
  }

  size_t Size() const { return pdes_.size(); }

 private:
  std::vector<std::unique_ptr<PDE>> pdes_;
  std::unordered_map<std::string, size_t> index_;
};

// Single standalone operation that runs all PDE source terms each step.
// Replaces CalciumRecovery + DermalO2Source as separate operations.
struct CompositeFieldOp : public StandaloneOperationImpl {
  BDM_OP_HEADER(CompositeFieldOp);

  explicit CompositeFieldOp(CompositeField* fields) : fields_(fields) {}

  void operator()() override {
    fields_->ApplyAllSources(Simulation::GetActive());
  }

 private:
  CompositeField* fields_;
};

}  // namespace skibidy
}  // namespace bdm

#endif  // FIELDS_COMPOSITE_FIELD_H_
