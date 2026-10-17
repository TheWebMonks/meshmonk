#ifndef CORRESPONDENCEFILTERFACTORY_HPP
#define CORRESPONDENCEFILTERFACTORY_HPP

#include "BaseCorrespondenceFilter.hpp"
#include "CorrespondenceFilter.hpp"
#include "SymmetricCorrespondenceFilter.hpp"
#include <cstddef>
#include <memory>

namespace registration {

// The single construction site for correspondence filters: symmetric
// (push/pull) or plain, configured with the shared parameters.
inline std::unique_ptr<BaseCorrespondenceFilter>
make_correspondence_filter(bool symmetric, size_t num_neighbours,
                           float flag_threshold, bool equalize_push_pull) {
  if (symmetric) {
    auto f = std::make_unique<SymmetricCorrespondenceFilter>();
    f->set_parameters(num_neighbours, flag_threshold, equalize_push_pull);
    return f;
  }
  auto f = std::make_unique<CorrespondenceFilter>();
  f->set_parameters(num_neighbours, flag_threshold);
  return f;
}

} // namespace registration

#endif // CORRESPONDENCEFILTERFACTORY_HPP
