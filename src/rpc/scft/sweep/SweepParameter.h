#ifndef RPC_SWEEP_PARAMETER_H
#define RPC_SWEEP_PARAMETER_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <iostream>
#include <string>
#include <pscf/sweep/ParameterType.h>

namespace Pscf {
namespace Rpc {

   template <int D> class System;

   /**
   * Class for storing data about an individual sweep parameter.
   *
   * This class stores the information required to sweep a single 
   * parameter value of any of several types.  The type of parameter
   * is indicated in the public interface and parameter file format
   * by a string identifier with any of several allowed values.
   * Each parameter is also identified by one or two associated index 
   * values, denoted here by id(0) and id(1), that specify the index 
   * or indices for a subobject or array element with which the 
   * parameter is associated applied. Allowed string representations 
   * and meanings of parameter types are given below, along with the 
   * meaning of any associated index value or pair of values.
   * To indicate the meaning of index values, we use mId to denote 
   * a monomer type index, pId to denote a polymer species index, 
   * bId to denote the index of a block within a polymer, sId to
   * denote a solvent species index, lId to denote a lattice parameter
   * index, and wId to denote a wall index:
   * \code
   *  | Type        | Meaning                            | id(0) | id(1)
   *  | ----------- | ---------------------------------- | ----- | -----
   *  | kuhn        | monomer segment length             | mId   |
   *  | chi         | Flory-Huggins parameter            | mId   | mId
   *  | block       | block length                       | pId   | bId
   *  | solvent     | solvent size                       | sId   |
   *  | phi_polymer | polymer volume fraction            | pId   |
   *  | mu_polymer  | polymer chemical potential         | pId   |
   *  | phi_solvent | solvent volume fraction            | sId   |
   *  | mu_solvent  | solvent chemical potential         | sId   |
   *  | cell_param  | lattice parameter                  | lId   |
   *  | chi_bottom  | thin film wall/monomer interaction | mId   | 
   *  | chi_top     | thin film wall/monomer interaction | mId   | 
   * \endcode
   * The two indices for a Flory-Huggins chi parameter refer to indices
   * in the chi matrix maintained by Interaction. Changes to element
   * chi(i, j) automatically also update chi(j, i) for i !\ j, thus
   * maintaining the symmetry of the matrix.
   *
   * Each SweepParameter also has a "change" value that gives the 
   * intended difference between the final and initial value of the 
   * parameter over the course of a sweep, corresponding to a change
   * sweep parameter s over the range [0,1]. The initial value of each 
   * parameter is obtained from a query of the state of the parent 
   * system at the beginning of a sweep, and thus does not need to 
   * be supplied as part of the text format for a SweepParameter.
   *
   * A SweepParameter<D> object is initialized by reading the parameter
   * type, index or index and change value from a parameter file as a
   * a single line.  An overloaded >> operator is defined that allows 
   * a SweepParameter<D> object named "parameter" to be read from an 
   * istream named "in" using the syntax "in >> parameter". 
   *
   * The text format for a parameter of a type that requires a single
   * index id(0) is:
   *
   *    type id(0) change 
   *
   * where type indicates a type string, id(0) is an integer index value,
   * and change is the a floating point value for the change in parameter
   * value. The corresponding format for a parameter that requires two 
   * indices (e.g., block or chi) is instead: "type id(0) id(1) change".
   *
   * Most of the parameter types in the table above are hard-coded into
   * this class, but certain classes can also add "specialized"
   * parameters to the list SweepTmpl::parameterTypes_ using the method
   * SweepTmpl::addParameterTypes. This list of specialized parameters 
   * is then consulted if the user-specified parameter type is not found
   * in the list of hard-coded parameter names. In order to add parameters
   * to SweepTmpl::parameterTypes_, a class must be a subclass of
   * ParameterModifier. 
   * 
   * Currently, the only classes that can add parameters to 
   * parameterTypes_ are Iterators. The Iterator for a given calculation
   * is determined via a Factory, so the type of Iterator is not known
   * until run time. Therefore, Iterators can add sweepable parameters
   * to parameterTypes_ as necessary at run time, depending on the type of
   * Iterator used in that calculation. As a design principle, an object
   * that is created by a Factory and used for SCFT calculations should
   * also be a subclass of ParameterModifier so that it can add sweepable
   * parameters at run time.
   *
   * \ingroup Rpc_Scft_Sweep_Module
   */
   template <int D>
   class SweepParameter
   {

   public:

      /**
      * Default constructor.
      */
      SweepParameter();

      /**
      * Constructor that stores a pointer to parent system.
      *
      * \param system  parent system
      */
      SweepParameter(System<D>& system);

      /**
      * Set the system associated with this object.
      *
      * Invoke this function on objects created with the default 
      * constructor to create an association with a parent system.
      *
      * \param system  parent system
      */
      void setSystem(System<D>& system)
      {  systemPtr_ = &system; }

      /**
      * Set the pointer to the array of specialized sweep parameter types.
      * 
      * \param array  array of specialized parameter types
      */
      void setParameterTypesArray(GArray<ParameterType>& array)
      {  parameterTypesPtr_ = &array; }

      /**
      * Get the ParameterType object for a specialized sweep parameter.
      * 
      * An error will result if this SweepParameter is not specialized
      * and this function is called.
      */
      ParameterType& parameterType() const;

      /**
      * Get the array index for the specialized sweep parameter.
      */
      int parameterTypeId() const
      {  return parameterTypeId_; }

      /**
      * Is this SweepParameter a specialized parameter type? 
      */
      bool isSpecialized() const
      {  return (parameterTypeId_ != -1); }

      /**
      * Get and store initial value of this parameter.
      *
      * This function is called before a sweep begins, and simply gets
      * current values of this parameter.
      */
      void getInitial();

      /**
      * Update the corresponding parameter value in the System.
      *
      * \param newVal  new value for this parameter (input)
      */
      void update(double newVal);

      /**
      * Return a string representation of the parameter type. 
      */
      std::string type() const;

      /**
      * Write the parameter type to an output stream.
      *
      * \param out  output file stream
      */
      void writeParamType(std::ostream& out) const;

      /**
      * Get id for a sub-object or element to which this is applied.
      *
      * This function returns a value from the id_ array. Elements of 
      * array store indices associating with a parameter value, which
      * may, for example, identify a molecule type index, one or more 
      * monomer type indices or a unit cell parameter index. 
      * 
      * Different types of parameters require either 1 or 2 such 
      * identifiers.  The number of required identifiers is returned 
      * by the function nId(). 
      *
      * \param i array index to access
      */
      int id(int i) const
      {
         UTIL_CHECK(i < nId_);  
         return id_[i];
      }

      /**
      * Number of indices associated with this type of parameter.
      *
      * See documentation of id(int).
      */
      int nId() const
      {  return nId_; }

      /**
      * Get the current system parameter value.
      */
      double current()
      {  return get_(); }

      /**
      * Get the initial system parameter value.
      */
      double initial() const
      {  return initial_; }

      /**
      * Get the total change planned for this parameter during sweep.
      */
      double change() const
      {  return change_; }

      /**
      * Serialize to or from an archive.
      *
      * \param ar Archive object
      * \param version archive format version index
      */
      template <class Archive>
      void serialize(Archive ar, const unsigned int version);

   private:

      /// Enumeration of allowed parameter types.
      enum ParamType { Block, Chi, Kuhn, Phi_Polymer, Phi_Solvent,
                       Mu_Polymer, Mu_Solvent, Solvent, Cell_Param, 
                       Special, Null};

      /// Type of parameter associated with an object of this class.
      ParamType type_;

      /// Number of identifiers needed for this parameter type.
      int nId_;

      /// Identifier indices (e.g., molecule or monomer type index)
      DArray<int> id_;

      /// Initial parameter value, retrieved from system at start of sweep.
      double initial_;

      /// Change in parameter
      double change_;

      /// Pointer to the parent system.
      System<D>* systemPtr_;

      /// Pointer to the parameterTypes_ array of the Sweep base class.
      GArray<ParameterType>* parameterTypesPtr_;

      /// Index of parameter type within parameterTypes_ array if found
      int parameterTypeId_;

      /**
      * Read type of parameter being swept, and set number of identifiers.
      *
      * \param in  input stream from param file.
      */
      void readParamType(std::istream& in);

      /**
      * Gets the current system parameter value.
      */
      double get_();

      /**
      * Set a new system parameter value.
      *
      * \param newVal  new value for this parameter.
      */
      void set_(double newVal);

   // friends:

      template <int U>
      friend
      std::istream& operator >> (std::istream&, SweepParameter<U>&);

      template <int U>
      friend
      std::ostream& operator << (std::ostream&, SweepParameter<U> const&);

   };

   /**
   * Inserter for reading a SweepParameter from an istream.
   *
   * \param in  input stream
   * \param param  SweepParameter<D> object to read
   */
   template <int D>
   std::istream& operator >> (std::istream& in, 
                              SweepParameter<D>& param);

   /**
   * Extractor for writing a SweepParameter to ostream.
   *
   * \param out  output stream
   * \param param  SweepParameter<D> object to write
   */
   template <int D>
   std::ostream& operator << (std::ostream& out, 
                              SweepParameter<D> const & param);

}
}
#include "SweepParameter.tpp"
#endif
