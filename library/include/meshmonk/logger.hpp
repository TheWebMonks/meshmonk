#ifndef MESHMONK_LOGGER_HPP
#define MESHMONK_LOGGER_HPP

#include <atomic>
#include <iostream>
#include <string>

namespace meshmonk {

enum class LogLevel {
  Silent,
  Error,
  Warning,
  Info,
  Debug,
};

namespace detail {

// Thread-safe atomic storage for the current log level
inline std::atomic<LogLevel> &global_log_level() {
  static std::atomic<LogLevel> level{LogLevel::Warning};
  return level;
}

} // namespace detail

inline void set_log_level(LogLevel level) {
  detail::global_log_level().store(level, std::memory_order_relaxed);
}

inline LogLevel get_log_level() {
  return detail::global_log_level().load(std::memory_order_relaxed);
}

// Log a message at the given level. No-op if level > get_log_level().
inline void log(LogLevel level, const std::string &msg) {
  if (level == LogLevel::Silent)
    return;
  if (static_cast<int>(level) > static_cast<int>(get_log_level()))
    return;

  const char *prefix = "";
  switch (level) {
  case LogLevel::Error:
    prefix = "[ERROR] ";
    break;
  case LogLevel::Warning:
    prefix = "[WARNING] ";
    break;
  case LogLevel::Info:
    prefix = "[INFO] ";
    break;
  case LogLevel::Debug:
    prefix = "[DEBUG] ";
    break;
  default:
    break;
  }
  std::cerr << prefix << msg << "\n";
}

} // namespace meshmonk

#endif // MESHMONK_LOGGER_HPP
