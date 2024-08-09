
template<class T>
template<class U> requires std::convertible_to<U,T>
LinearCombination<T>::LinearCombination(U value) : bias(value) {}

template<class T>
LinearCombination<T>::LinearCombination(std::initializer_list<Term> const& lst) : terms(lst.size())
{
    std::copy(lst.begin(), lst.end(), terms.begin());
}


template<class T>
LinearCombination<T> operator+(LinearCombination<T> lhs, LinearCombination<T> const& rhs)
{
    return lhs += rhs;
}


template<class T>
LinearCombination<T> operator-(LinearCombination<T> lhs, LinearCombination<T> const& rhs)
{
    return lhs -= rhs;
}

template<class T>
LinearCombination<T> operator-(LinearCombination<T> lc)
{
    return lc *= -1;
}


template<class T>
LinearCombination<T> operator*(LinearCombination<T> lhs, Scalar rhs)
{
    return lhs *= rhs;
}


template<class T>
LinearCombination<T> operator*(Scalar lhs, LinearCombination<T> rhs)
{
    return rhs *= lhs;
}


template<class T>
LinearCombination<T>& operator*=(LinearCombination<T>& lhs, Scalar rhs)
{
    if (rhs == 0)
        return lhs = LinearCombination<T>();

    for (Term& term : lhs.terms)
        term.coeff *= rhs;
    lhs.bias *= rhs;

    return lhs;
}

template<class T>
LinearCombination<T>& LinearCombination<T>::operator-=(LinearCombination<T> const& rhs)
{
    *this *= -1;
    *this += rhs;
    return *this *= -1;
}


template<class T>
LinearCombination<T>& LinearCombination<T>::operator+=(LinearCombination<T> const& rhs)
{
    terms.reserve(terms.size() + rhs.terms.size());

    for (Term term : rhs.terms)
    {
        auto it = std::find(terms.begin(), terms.end(), term);
        if (it != terms.end())
            it->coeff += term.coeff;
        else
            terms.push_back(term);
    }
    bias += rhs.bias;

    return *this;
}

template<class T> 
template <class Field>
T LinearCombination<T>::evaluate(Field const& field) const
requires requires(Field field, Index idx)
{
    getFieldValue(field, idx);
    {getFieldValue(field, idx)} -> std::same_as<T>;
}
{
    T result = bias;
    for (auto [coeff, idx] : terms)
        result += coeff * getFieldValue(field, idx);
    return result;
}

template<class T>
std::ostream& operator<<(std::ostream& os, LinearCombination<T> const& lc)
{
    if (lc.terms.size() == 0)
        return os << lc.bias << '\n';
    
    if (lc.terms.front().coeff < 0)
        os << '-';

    auto abs = [](Scalar c) {return (c > 0 ? c : -c);};
    for (Index i = 0; i < lc.terms.size()-1; i++)
        os << abs(lc.terms[i].coeff) << "*a" << lc.terms[i].idx << (lc.terms[i+1].coeff > 0 ? " + " : " - ");
    os << abs(lc.terms.back().coeff) << "*a" << lc.terms.back().idx;

    if (lc.bias != T())
        os << " + " << lc.bias;
    
    return os << '\n';
}