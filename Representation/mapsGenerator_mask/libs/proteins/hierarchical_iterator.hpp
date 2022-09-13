#pragma once


#include <vector>
#include <functional>


/* The First Level Iterator */
template <class T>
class FirstLevelIteratorGenerator {
  class parent_iterator;

 public:
  explicit FirstLevelIteratorGenerator(const std::vector<T*> &vector)
      : vector_(vector) {}

  class iterator;
  iterator begin() { return iterator(vector_.begin()); }
  iterator end() { return iterator(vector_.end()); }

  class const_iterator;
  const_iterator begin() const { return const_iterator(vector_.begin()); }
  const_iterator end() const { return const_iterator(vector_.end()); }

 private:
  const std::vector<T*> &vector_;
};

template <class T>
class FirstLevelIteratorGenerator<T>::parent_iterator {
 public:
  parent_iterator() {}
  explicit parent_iterator(const typename std::vector<T*>::const_iterator &it) : it_(it) {}
  virtual ~parent_iterator() {}
  parent_iterator& operator=(const parent_iterator &other) = default;
  bool operator==(const parent_iterator &other) const  { return it_ == other.it_; }
  bool operator!=(const parent_iterator &other) const { return it_ != other.it_; }
  parent_iterator& operator++() { ++it_; return *this; }

 protected:
  typename std::vector<T*>::const_iterator it_;
};

template <class T>
class FirstLevelIteratorGenerator<T>::iterator : public parent_iterator {
 public:
  iterator() {}
  explicit iterator(const typename std::vector<T*>::const_iterator &it)
      : parent_iterator(it) {}
  T& operator*() const { return **this->it_; }
  T* operator->() const { return *this->it_; }
};

template <class T>
class FirstLevelIteratorGenerator<T>::const_iterator : public parent_iterator {
 public:
  const_iterator() {}
  explicit const_iterator(const typename std::vector<T*>::const_iterator &it)
      : parent_iterator(it) {}
  const T& operator*() const { return **this->it_; }
  const T* operator->() const { return *this->it_; }
};


/* High Level Iterator */
template <
  class T,
  class G,
  class ItGen = FirstLevelIteratorGenerator<G>
>
class HighLevelIteratorGenerator {
  class parent_iterator;
  typedef std::function<ItGen(const T&)> getIterGenType;

 public:
  HighLevelIteratorGenerator(const std::vector<T*> &vector,
                             const getIterGenType &get_first_level_iter_gen)
      : vector_(vector), get_first_level_iter_gen_(get_first_level_iter_gen) {}

  class iterator;
  iterator begin() { return iterator(vector_, get_first_level_iter_gen_, false); }
  iterator end() { return iterator(vector_, get_first_level_iter_gen_, true); }

  class const_iterator;
  const_iterator begin() const { return const_iterator(vector_, get_first_level_iter_gen_, false); }
  const_iterator end() const { return const_iterator(vector_, get_first_level_iter_gen_, true); }

 private:
  const std::vector<T*> &vector_;
  getIterGenType get_first_level_iter_gen_;
};


template <class T, class G, class ItGen>
class HighLevelIteratorGenerator<T, G, ItGen>::parent_iterator {
 public:
  parent_iterator() {}
  parent_iterator(const std::vector<T*> &vector,
                  const getIterGenType &get_first_level_iter_gen,
                  bool end_iterator);
  virtual ~parent_iterator() {}
  parent_iterator& operator=(const parent_iterator &it) = default;
  bool operator==(const parent_iterator &it) const;
  bool operator!=(const parent_iterator &it) const;
  parent_iterator& operator++();

 protected:
  size_t i_;
  const std::vector<T*> *vector_;
  typename ItGen::iterator first_level_it_;
  getIterGenType get_first_level_iter_gen_;
};

template <class T, class G, class ItGen>
HighLevelIteratorGenerator<T, G, ItGen>::parent_iterator::parent_iterator(
                          const std::vector<T*> &vector,
                          const getIterGenType &get_first_level_iter_gen,
                          bool end_iterator)
    : i_(0), vector_(&vector), get_first_level_iter_gen_(get_first_level_iter_gen) {
  if (end_iterator || !vector_->size()) {
    i_ = vector_->size();
    return;
  }
  while (true) {
    while (i_ < vector_->size() && !vector_->at(i_))
      i_++;
    if (i_ == vector_->size())
      return;

    first_level_it_ = get_first_level_iter_gen_(*vector_->at(i_)).begin();
    auto first_level_it_end = get_first_level_iter_gen_(*vector_->at(i_)).end();
    if (first_level_it_ != first_level_it_end) {
      return;
    } else {
      i_++;
    }
  }
}

template <class T, class G, class ItGen>
bool HighLevelIteratorGenerator<T, G, ItGen>::parent_iterator::operator==(
        const parent_iterator &other) const {
  return i_ == other.i_ && (i_ == vector_->size() || first_level_it_ == other.first_level_it_);
}

template <class T, class G, class ItGen>
bool HighLevelIteratorGenerator<T, G, ItGen>::parent_iterator::operator!=(
        const parent_iterator &other) const {
  return i_ != other.i_ || (i_ < vector_->size() && first_level_it_ != other.first_level_it_);
}

template <class T, class G, class ItGen>
typename HighLevelIteratorGenerator<T, G, ItGen>::parent_iterator&
 HighLevelIteratorGenerator<T, G, ItGen>::parent_iterator::operator++() {
  if (++first_level_it_ != get_first_level_iter_gen_(*vector_->at(i_)).end())
    return *this;
  do {
    do {
      i_++;
    } while (i_ < vector_->size() && !vector_->at(i_));
    if (i_ == vector_->size()) {
      return *this;
    }
    first_level_it_ = get_first_level_iter_gen_(*vector_->at(i_)).begin();
  } while (first_level_it_ == get_first_level_iter_gen_(*vector_->at(i_)).end());
  return *this;
}


template <class T, class G, class ItGen>
class HighLevelIteratorGenerator<T, G, ItGen>::iterator : public parent_iterator {
 public:
  iterator() {}
  iterator(const std::vector<T*> &vector,
           const getIterGenType &get_first_level_iter_gen,
           bool end_iterator)
      : parent_iterator(vector, get_first_level_iter_gen, end_iterator) {}
  G& operator*() const { return *this->first_level_it_; }
  G* operator->() const { return &(*this->first_level_it_); }
};


template <class T, class G, class ItGen>
class HighLevelIteratorGenerator<T, G, ItGen>::const_iterator : public parent_iterator {
 public:
  const_iterator() {}
  const_iterator(const std::vector<T*> &vector,
                 const getIterGenType &get_first_level_iter_gen,
                 bool end_iterator)
      : parent_iterator(vector, get_first_level_iter_gen, end_iterator) {}
  const G& operator*() const { return *this->first_level_it_; }
  const G* operator->() const { return &(*this->first_level_it_); }
};
