export Dimension, Degree, MultiIndex, Form

const Dimension = Integer
const Degree = Integer
const MultiIndex{D} = NTuple{D,Integer}
const Form = Integer

abstract type AbstractPatch end
abstract type AbstractBox{D} <: AbstractPatch end
abstract type AbstractSimplex{D} <: AbstractPatch end
abstract type AbstractPrism{D1,D2} <: AbstractPatch end
