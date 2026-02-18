#ifndef CHECKPOINT_H_
#define CHECKPOINT_H_

// Lightweight checkpoint: save/load DiffusionGrid data and agent state
// to a directory of binary files. Used for fork-based batch runs where
// configs share a common untreated prefix up to some day.
//
// Env vars:
//   SKIBIDY_CKPT_SAVE_DIR  - directory to save checkpoint into
//   SKIBIDY_CKPT_LOAD_DIR  - directory to load checkpoint from
//   SKIBIDY_CKPT_STEP      - step number for save (run to this step, save, exit)
//                             For load: the step is read from the checkpoint.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <sys/stat.h>

#include "biodynamo.h"
#include "core/field_names.h"

namespace bdm {
namespace skibidy {
namespace checkpoint {

// Header written to each grid file for validation
struct GridFileHeader {
  char magic[8] = {'S', 'K', 'C', 'K', 'P', 'T', '0', '1'};
  uint64_t num_boxes = 0;
  uint64_t resolution = 0;
  int32_t dims[6] = {};
};

// Header for agent state file
struct AgentFileHeader {
  char magic[8] = {'S', 'K', 'A', 'G', 'N', 'T', '0', '1'};
  uint64_t num_agents = 0;
};

// Per-agent record (position + key state)
struct AgentRecord {
  double x, y, z;
  double diameter;
  int32_t type;      // 0=keratinocyte, 1=immune, 2=tumor, 3=fibroblast
  int32_t subtype;   // stratum for kerato, immune_type for immune, etc.
  double age;
};

inline bool MkdirP(const std::string& path) {
  return mkdir(path.c_str(), 0755) == 0 || errno == EEXIST;
}

/// Save all DiffusionGrid data to binary files in the given directory.
/// Also saves step count and agent state.
inline bool SaveCheckpoint(Simulation* sim, uint64_t step,
                           const std::string& dir) {
  if (!MkdirP(dir)) {
    Log::Error("Checkpoint", "Cannot create directory: ", dir);
    return false;
  }

  auto* rm = sim->GetResourceManager();

  // Save step number
  {
    std::string path = dir + "/step.bin";
    std::ofstream f(path, std::ios::binary);
    f.write(reinterpret_cast<const char*>(&step), sizeof(step));
  }

  // Save each DiffusionGrid
  rm->ForEachDiffusionGrid([&](DiffusionGrid* grid) {
    std::string name = grid->GetSubstanceName();
    std::string path = dir + "/grid_" + name + ".bin";
    std::ofstream f(path, std::ios::binary);

    GridFileHeader hdr;
    hdr.num_boxes = grid->GetNumBoxes();
    hdr.resolution = grid->GetResolution();
    auto dims = grid->GetDimensions();
    for (int i = 0; i < 6; i++) hdr.dims[i] = dims[i];

    f.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));

    const real_t* data = grid->GetAllConcentrations();
    f.write(reinterpret_cast<const char*>(data),
            hdr.num_boxes * sizeof(real_t));
  });

  // Save agent count and basic state
  {
    std::string path = dir + "/agents.bin";
    std::ofstream f(path, std::ios::binary);

    // Count agents first
    uint64_t count = 0;
    rm->ForEachAgent([&](Agent*) { count++; });

    AgentFileHeader hdr;
    hdr.num_agents = count;
    f.write(reinterpret_cast<const char*>(&hdr), sizeof(hdr));

    // Write each agent's position and diameter (enough for re-creation)
    rm->ForEachAgent([&](Agent* agent) {
      AgentRecord rec;
      auto pos = agent->GetPosition();
      rec.x = pos[0];
      rec.y = pos[1];
      rec.z = pos[2];
      rec.diameter = agent->GetDiameter();
      rec.type = 0;
      rec.subtype = 0;
      rec.age = 0;
      f.write(reinterpret_cast<const char*>(&rec), sizeof(rec));
    });
  }

  Log::Info("Checkpoint", "Saved checkpoint at step ", step, " to ", dir,
            " (", rm->GetNumAgents(), " agents)");
  return true;
}

/// Load DiffusionGrid data from checkpoint directory.
/// Overwrites current grid concentrations. Returns the saved step number.
inline uint64_t LoadCheckpoint(Simulation* sim, const std::string& dir) {
  auto* rm = sim->GetResourceManager();

  // Load step number
  uint64_t step = 0;
  {
    std::string path = dir + "/step.bin";
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) {
      Log::Fatal("Checkpoint", "Cannot open step file: ", path);
    }
    f.read(reinterpret_cast<char*>(&step), sizeof(step));
  }

  // Load each DiffusionGrid
  rm->ForEachDiffusionGrid([&](DiffusionGrid* grid) {
    std::string name = grid->GetSubstanceName();
    std::string path = dir + "/grid_" + name + ".bin";
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) {
      Log::Warning("Checkpoint", "Grid file not found: ", path,
                   " (skipping, grid stays at initial state)");
      return;
    }

    GridFileHeader hdr;
    f.read(reinterpret_cast<char*>(&hdr), sizeof(hdr));

    // Validate magic
    if (std::memcmp(hdr.magic, "SKCKPT01", 8) != 0) {
      Log::Fatal("Checkpoint", "Bad magic in ", path);
    }

    // Validate grid size matches
    if (hdr.num_boxes != grid->GetNumBoxes()) {
      Log::Warning("Checkpoint", "Grid size mismatch for ", name,
                   ": checkpoint has ", hdr.num_boxes, " boxes, current has ",
                   grid->GetNumBoxes(), " (skipping)");
      return;
    }

    // Read concentrations into the grid.
    // We zero out first, then add saved values via ChangeConcentrationBy.
    // This respects lower/upper thresholds set during field registration.
    const real_t* current = grid->GetAllConcentrations();
    for (size_t i = 0; i < hdr.num_boxes; i++) {
      real_t saved_val;
      f.read(reinterpret_cast<char*>(&saved_val), sizeof(saved_val));
      real_t delta = saved_val - current[i];
      if (std::abs(delta) > 1e-15) {
        grid->ChangeConcentrationBy(i, delta);
      }
    }
  });

  Log::Info("Checkpoint", "Loaded checkpoint from ", dir, " at step ", step);
  return step;
}

/// Check env vars and return checkpoint config.
struct CheckpointConfig {
  bool save = false;
  bool load = false;
  std::string save_dir;
  std::string load_dir;
  uint64_t save_step = 0;
};

inline CheckpointConfig GetCheckpointConfig() {
  CheckpointConfig cfg;

  const char* save_dir = std::getenv("SKIBIDY_CKPT_SAVE_DIR");
  const char* load_dir = std::getenv("SKIBIDY_CKPT_LOAD_DIR");
  const char* step_str = std::getenv("SKIBIDY_CKPT_STEP");

  if (save_dir && save_dir[0]) {
    cfg.save = true;
    cfg.save_dir = save_dir;
    cfg.save_step = step_str ? std::atoi(step_str) : 0;
  }
  if (load_dir && load_dir[0]) {
    cfg.load = true;
    cfg.load_dir = load_dir;
  }

  return cfg;
}

}  // namespace checkpoint
}  // namespace skibidy
}  // namespace bdm

#endif  // CHECKPOINT_H_
