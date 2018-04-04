/*
Copyright (c) 2013, Ryan Haining, Aaron Josephs, Google
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

  Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

  Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef ITER_PRODUCT_HPP_
#define ITER_PRODUCT_HPP_

#include "internal/iterbase.hpp"

#include <iterator>
#include <tuple>
#include <utility>
#include <array>

namespace iter {
  namespace impl {
    template <typename... Containers>
    class Productor;

    template <typename Container, typename... RestContainers>
    class Productor<Container, RestContainers...>;

    template <>
    class Productor<>;
  }

  template <typename... Containers>
  impl::Productor<Containers...> product(Containers&&...);
}

// specialization for at least 1 template argument
template <typename Container, typename... RestContainers>
class iter::impl::Productor<Container, RestContainers...> {
  friend Productor iter::product<Container, RestContainers...>(
      Container&&, RestContainers&&...);

  template <typename... RC>
  friend class Productor;

  using ProdIterDeref =
      std::tuple<iterator_deref<Container>, iterator_deref<RestContainers>...>;

 private:
  Container container;
  Productor<RestContainers...> rest_products;
  Productor(Container&& in_container, RestContainers&&... rest)
      : container(std::forward<Container>(in_container)),
        rest_products{std::forward<RestContainers>(rest)...} {}

 public:
  Productor(Productor&&) = default;
  class Iterator
      : public std::iterator<std::input_iterator_tag, ProdIterDeref> {
   private:
    using RestIter = typename Productor<RestContainers...>::Iterator;

    iterator_type<Container> iter;
    iterator_type<Container> begin;

    RestIter rest_iter;
    RestIter rest_end;

   public:
    constexpr static const bool is_base_iter = false;
    Iterator(const iterator_type<Container>& it, RestIter&& rest,
        RestIter&& in_rest_end)
        : iter{it}, begin{it}, rest_iter{rest}, rest_end{in_rest_end} {}

    void reset() {
      this->iter = this->begin;
    }

    Iterator& operator++() {
      ++this->rest_iter;
      if (!(this->rest_iter != this->rest_end)) {
        this->rest_iter.reset();
        ++this->iter;
      }
      return *this;
    }

    Iterator operator++(int) {
      auto ret = *this;
      ++*this;
      return ret;
    }

    bool operator!=(const Iterator& other) const {
      return this->iter != other.iter
             && (RestIter::is_base_iter || this->rest_iter != other.rest_iter);
    }

    bool operator==(const Iterator& other) const {
      return !(*this != other);
    }

    ProdIterDeref operator*() {
      return std::tuple_cat(
          std::tuple<iterator_deref<Container>>{*this->iter}, *this->rest_iter);
    }

    ArrowProxy<ProdIterDeref> operator->() {
      return {**this};
    }
  };

  Iterator begin() {
    return {std::begin(this->container), std::begin(this->rest_products),
        std::end(this->rest_products)};
  }

  Iterator end() {
    return {std::end(this->container), std::end(this->rest_products),
        std::end(this->rest_products)};
  }
};

template <>
class iter::impl::Productor<> {
 public:
  Productor(Productor&&) = default;
  class Iterator : public std::iterator<std::input_iterator_tag, std::tuple<>> {
   public:
    constexpr static const bool is_base_iter = true;

    void reset() {}

    Iterator& operator++() {
      return *this;
    }

    Iterator operator++(int) {
      auto ret = *this;
      ++*this;
      return ret;
    }

    // see note in zip about base case operator!=
    bool operator!=(const Iterator&) const {
      return false;
    }

    bool operator==(const Iterator& other) const {
      return !(*this != other);
    }

    std::tuple<> operator*() const {
      return {};
    }
  };

  Iterator begin() {
    return {};
  }

  Iterator end() {
    return {};
  }
};

template <typename... Containers>
iter::impl::Productor<Containers...> iter::product(Containers&&... containers) {
  return {std::forward<Containers>(containers)...};
}

namespace iter {
  constexpr std::array<std::tuple<>, 1> product() {
    return {{}};
  }
}

#endif
