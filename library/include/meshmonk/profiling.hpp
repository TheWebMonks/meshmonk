#pragma once

// Unconditionally include types used in both ifdef and else branches
#include <cstdint>
#include <string>
#include <unordered_map>

#ifdef MESHMONK_PROFILING

#include <chrono>

// Forward declaration — allows Profiler to declare scoped() before ScopedTimer
// is fully defined.
class ScopedTimer;

struct ProfileEntry {
  uint64_t total_us{0};
  uint64_t count{0};
};

class Profiler {
public:
  void record(const std::string& label, uint64_t us) {
    auto& e = entries_[label];
    e.total_us += us;
    e.count += 1;
  }
  // Out-of-line definition below, after ScopedTimer is fully defined.
  ScopedTimer scoped(const std::string& label);

  std::unordered_map<std::string, ProfileEntry> snapshot() const {
    return entries_;
  }
  void reset() { entries_.clear(); }
  bool empty() const { return entries_.empty(); }

  // Out-of-line definition below, after ScopedTimer is fully defined.
  uint64_t calibrate(size_t n);

private:
  std::unordered_map<std::string, ProfileEntry> entries_;
};

class ScopedTimer {
public:
  explicit ScopedTimer(const std::string& label, Profiler& p)
      : label_(label), profiler_(p),
        start_(std::chrono::steady_clock::now()) {}
  ~ScopedTimer() {
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::steady_clock::now() - start_).count();
    profiler_.record(label_, static_cast<uint64_t>(elapsed));
  }
  ScopedTimer(const ScopedTimer&) = delete;
  ScopedTimer& operator=(const ScopedTimer&) = delete;
  ScopedTimer(ScopedTimer&&) = delete;
  ScopedTimer& operator=(ScopedTimer&&) = delete;

private:
  std::string label_;
  Profiler& profiler_;
  std::chrono::steady_clock::time_point start_;
};

// Out-of-line definitions of Profiler methods that use ScopedTimer — required
// because ScopedTimer was only forward-declared when Profiler was defined above.
inline ScopedTimer Profiler::scoped(const std::string& label) {
  return ScopedTimer(label, *this);
}

// Run n empty ScopedTimer scopes and return nanoseconds-per-scope.
// Uses wall time of the full loop divided by n.
// The '__calibration__' entries are erased from the accumulator before
// returning so the calibration loop does not pollute profile data.
inline uint64_t Profiler::calibrate(size_t n) {
  auto t0 = std::chrono::steady_clock::now();
  for (size_t i = 0; i < n; ++i) {
    ScopedTimer t("__calibration__", *this);
  }
  auto t1 = std::chrono::steady_clock::now();
  entries_.erase("__calibration__");
  auto total_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
      t1 - t0).count();
  return static_cast<uint64_t>(
      total_ns / static_cast<int64_t>(n));
}

inline Profiler g_profiler;

#else  // MESHMONK_PROFILING not defined — no-op stubs

// These stubs ensure callsite patterns `auto _t = g_profiler.scoped("label");`
// and `g_profiler.calibrate(N)` compile with zero overhead.
struct ProfileEntry {
  uint64_t total_us{0};
  uint64_t count{0};
};
class ScopedTimer {};
class Profiler {
public:
  ScopedTimer scoped(const std::string&) { return {}; }
  void reset() {}
  bool empty() const { return true; }
  uint64_t calibrate(size_t) { return 0; }
  std::unordered_map<std::string, ProfileEntry> snapshot() const { return {}; }
};
inline Profiler g_profiler;

#endif  // MESHMONK_PROFILING
