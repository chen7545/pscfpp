#ifndef PSPC_SYSTEM_H
#define PSPC_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/


#include <util/param/ParamComposite.h>     // base class

#include <pspc/mcmove/McSimulator.h>     // member

#include <pspc/solvers/Mixture.h>          // member
#include <pspc/field/Domain.h>             // member
#include <pspc/field/FieldIo.h>            // member
#include <pspc/field/WFieldContainer.h>    // member
#include <pspc/field/CFieldContainer.h>    // member
#include <pspc/field/Mask.h>               // member
#include <pspc/field/RField.h>             // member
#include <pspc/field/RFieldDft.h>          // member

#include <pscf/homogeneous/Mixture.h>      // member

#include <util/random/Random.h>            // member
#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // member template
#include <util/containers/FSArray.h>       // member template

namespace Pscf {

   class Interaction;

namespace Pspc
{

   template <int D> class Iterator;
   template <int D> class IteratorFactory;
   template <int D> class Sweep;
   template <int D> class SweepFactory;
   template <int D> class Compressor;
   template <int D> class CompressorFactory;
   template <int D> class McMove;
   template <int D> class McMoveFactory;

   using namespace Util;

   /**
   * Main class for SCFT simulation of one system.
   *
   * A System has (among other components):
   *
   *    - a Mixture (a container for polymer and solvent solvers)
   *    - a Domain (a description of the crystal domain and discretization)
   *    - an Iterator
   *    - Monomer chemical fields in both basis and grid formats
   *    - Monomer concentration fields in both basis and grid formats
   *
   * In a parameter file format, the main block is a System{...} block that
   * contains subblocks for sub-objects.
   *
   * \ref user_param_pc_page "Parameter File Format"
   * \ingroup Pscf_Pspc_Module
   */
   template <int D>
   class System : public ParamComposite
   {

   public:

      /// Field defined on an real-space grid (an r-grid)
      typedef RField<D> Field;

      /// \name Construction and Destruction
      //@{

      /**
      * Constructor.
      */
      System();

      /**
      * Destructor.
      */
      ~System();

      //@}
      /// \name Lifetime (Actions)
      //@{

      /**
      * Process command line options.
      */
      void setOptions(int argc, char **argv);

      /**
      * Read input parameters (with opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);

      /**
      * Read input parameters from default param file.
      */
      void readParam();

      /**
      * Read body of parameter block (without opening and closing lines).
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Read command script from a file.
      *
      * \param in command script file.
      */
      void readCommands(std::istream& in);

      /**
      * Read commands from default command file.
      */
      void readCommands();

      //@}
      /// \name W Field Modifiers
      //@{

      /**
      * Read chemical potential fields in symmetry adapted basis format.
      *
      * This function opens and reads the file with the name given by the
      * "filename" string, which must contain chemical potential fields
      * in symmetry-adapted basis format. The function copies these
      * fields to set new values for the system w fields in basis format, 
      * and also computes and resets the system w fields in r-space 
      * format. On exit, both w().basis() and w().rgrid() have been 
      * reset, w().hasData and w().isSymmetric() are true, and 
      * hasCFields() is false.
      *
      * \param filename name of input w-field basis file
      */
      void readWBasis(const std::string & filename);

      /**
      * Read chemical potential fields in real space grid (r-grid) format.
      *
      * This function opens and reads the file with the name given by the
      * "filename" string, which must contains chemical potential fields
      * in real space grid (r-grid) format. The function sets values for 
      * system w fields in r-grid format. It does not set attempt to set
      * field values in symmetry-adapted basis format, because it cannot
      * be known whether the r-grid field exhibits the declared space
      * group symmetry.  On exit, w().rgrid() is reset and w().hasData() 
      * is true, while w().isSymmetric() and hasCFields() are false.
      *
      * \param filename  name of input w-field basis file
      */
      void readWRGrid(const std::string & filename);

      /**
      * Set chemical potential fields, in symmetry-adapted basis format.
      *
      * This function sets values for w fields in both symmetry adapted 
      * and r-grid format.  On exit, values of both w().basis() and 
      * w().rgrid() are reset, w().hasData() and w().isSymmetric() are 
      * true, and hasCFields() is false. 
      *
      * \param fields  array of new w (chemical potential) fields
      */
      void setWBasis(DArray< DArray<double> > const & fields);

      /**
      * Set new w fields, in real-space (r-grid) format.
      *
      * This function set values for w fields in r-grid format, but does
      * not set components the symmetry-adapted basis format.  On exit, 
      * w.rgrid() is reset, w().hasData() is true, hasCFields() is false 
      * and w().isSymmetric() is false.
      *
      * \param fields  array of new w (chemical potential) fields
      */
      void setWRGrid(DArray< RField<D> > const & fields);

      //@}
      /// \name Unit Cell Modifiers
      //@{

      /**
      * Set parameters of the associated unit cell.
      *
      * \param unitCell  new UnitCell<D> (i.e., new parameters)
      */
      void setUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set state of the associated unit cell.
      *
      * \param lattice  lattice system
      * \param parameters  array of new unit cell parameters.
      */
      void setUnitCell(typename UnitCell<D>::LatticeSystem lattice,
                       FSArray<double, 6> const & parameters);

      /**
      * Set parameters of the associated unit cell.
      *
      * \param parameters  array of new unit cell parameters.
      */
      void setUnitCell(FSArray<double, 6> const & parameters);

      //@}
      /// \name Primary Field Theory Computations
      //@{

      /**
      * Solve the modified diffusion equation once, without iteration.
      *
      * This function calls the Mixture::compute() function to solve
      * the statistical mechanics problem for a non-interacting system
      * subjected to the currrent chemical potential fields. This
      * requires solution of the modified diffusion equation for all 
      * polymers, computation of Boltzmann weights for all solvents, 
      * computation of molecular partition * functions for all species, 
      * and computation of concentration fields for blocks and solvents, 
      * and computation of overall concentrations for all monomer types. 
      * This function does not compute the canonical (Helmholtz) free 
      * energy or grand-canonical free energy (i.e., pressure). Upon 
      * return, the flag hasCFields is set true.
      *
      * If argument needStress == true, then this function also calls
      * Mixture<D>::computeStress() to compute the stress.
      *
      * \param needStress true if stress is needed, false otherwise
      */
      void compute(bool needStress = false);

      /**
      * Iteratively solve a SCFT problem.
      *
      * This function calls the iterator to attempt to solve the SCFT
      * problem for the current mixture and system parameters, using
      * the current chemical potential fields and current unit cell 
      * parameter values as initial guesses.  On exit, hasCFields is 
      * set true whether or not convergence is obtained to within the 
      * desired tolerance.  The Helmholtz free energy and pressure are 
      * computed if and only if convergence is obtained.
      *
      * \pre The w().hasData() flag must be true on entry.
      *
      * \param isContinuation true if continuation within a sweep.
      * \return returns 0 for successful convergence, 1 for failure.
      */
      int iterate(bool isContinuation = false);

      /**
      * Sweep in parameter space, solving an SCF problem at each point.
      *
      * This function uses a Sweep object that was initialized in the
      * parameter file to solve the SCF problem at a sequence of points
      * along a line in parameter space. The nature of this sequence of
      * points is determined by implementation of a subclass of Sweep
      * and the parameters passed to the sweep object in the parameter
      * file.  The Iterator that is initialized in the parameter file
      * is called at each state point.
      *
      * An Exception is thrown if this is called when no Sweep has been
      * created (i.e., if hasSweep() == false).
      */
      void sweep();

      /**
      * Perform a field theoretic Monte-Carlo simulation.
      *
      * Perform a field theoretic Monte-Carlo simulation using the 
      * partial saddle-point approximation. 
      * 
      * \param nStep  number of Monte-Carlo steps
      */
      void simulate(int nStep);

      //@}
      /// \name Thermodynamic Properties
      //@{

      /**
      * Compute free energy density and pressure for current fields.
      *
      * This function should be called after a successful call of
      * Iterator::solve(). Resulting values are stored and then
      * accessed by the fHelmholtz() and pressure() functions.
      */
      void computeFreeEnergy();

      /**
      * Get precomputed Helmoltz free energy per monomer / kT.
      *
      * The value retrieved by this function is computed by the
      * computeFreeEnergy() function.
      */
      double fHelmholtz() const;

      /**
      * Get precomputed pressure x monomer volume kT.
      *
      * The value retrieved by this function is computed by the
      * computeFreeEnergy() function.
      */
      double pressure() const;

      //@}
      /// \name Output Operations (correspond to command file commands)
      //@{

      /**
      * Write chemical potential fields in symmetrized basis format.
      *
      * \param filename name of output file
      */
      void writeWBasis(const std::string & filename) const;

      /**
      * Write chemical potential fields in real space grid (r-grid) format.
      *
      * \param filename name of output file
      */
      void writeWRGrid(const std::string & filename) const;

      /**
      * Write concentration fields in symmetrized basis format.
      *
      * \param filename name of output file
      */
      void writeCBasis(const std::string & filename) const;

      /**
      * Write concentration fields in real space grid (r-grid) format.
      *
      * \param filename name of output file
      */
      void writeCRGrid(const std::string & filename) const;

      /**
      * Write c-fields for all blocks and solvents in r-grid format.
      *
      * Writes concentrations for all blocks of all polymers and all
      * solvent species in r-grid format. Columns associated with blocks
      * appear ordered by polymer id and then by block id, followed by
      * solvent species ordered by solvent id.
      *
      * \param filename name of output file
      */
      void writeBlockCRGrid(const std::string & filename) const;

      /**
      * Write slice of a propagator at fixed s in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      * \param segmentId  integer integration step index
      */
      void writeQSlice(std::string const & filename,
                       int polymerId, int blockId,
                       int directionId, int segmentId)  const;

      /**
      * Write the final slice of a propagator in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      */
      void writeQTail(std::string const & filename, int polymerId,
                      int blockId, int directionId)  const;

      /**
      * Write one propagator for one block, in r-grid format.
      *
      * \param filename  name of output file
      * \param polymerId  integer id of the polymer
      * \param blockId  integer id of the block within the polymer
      * \param directionId  integer id of the direction (0 or 1)
      */
      void writeQ(std::string const & filename, int polymerId,
                  int blockId, int directionId)  const;

      /**
      * Write all propagators of all blocks, each to a separate file.
      *
      * Write all propagators for both directions for all blocks
      * of all polymers, with each propagator in a separate file.
      * The function writeQ is called internally for each propagator,
      * and is passed an automatically generated file name.  The file
      * name for each propagator is given by a string of the form
      * (basename)_(ip)_(ib)_(id), where (basename) denotes the value
      * of the std::string function parameter basename, and where
      * (ip), (ib), and (id) denote the string representations of
      * a polymer indiex ip, a block index ib, and direction index id,
      * with id = 0 or 1. For example, if basename == "out/q", then
      * the file name of the propagator for direction 1 of block 2
      * of polymer 0 would be "out/q_0_2_1".
      *
      * \param basename  common prefix for output file names
      */
      void writeQAll(std::string const & basename);

      /**
      * Write parameter file to an ostream, omitting any sweep block.
      *
      * \param out output stream
      */
      void writeParamNoSweep(std::ostream& out) const;

      /**
      * Write thermodynamic properties to a file.
      *
      * This function outputs Helmholtz free energy per monomer,
      * pressure (in units of kT per monomer volume), and the
      * volume fraction and chemical potential of each species.
      *
      * \param out output stream
      */
      void writeThermo(std::ostream& out) const;

      /**
      * Output information about stars and symmetrized basis functions.
      *
      * This function opens a file with the specified filename and then
      * calls Basis::outputStars.
      *
      * \param outFileName name of output file
      */
      void writeStars(const std::string & outFileName) const;

      /**
      * Output information about waves.
      *
      * This function opens a file with the specified filename and then
      * calls Basis::outputWaves.
      *
      * \param outFileName name of output file for wave data
      */
      void writeWaves(const std::string & outFileName) const;

      //@}
      /// \name Field Operations (correspond to command file commands)
      //@{

      /**
      * Convert a field from symmetrized basis format to r-grid format.
      *
      * \param inFileName name of input file (basis format)
      * \param outFileName name of output file (r-grid format)
      */
      void basisToRGrid(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Convert a field from real-space grid to symmetrized basis format.
      *
      * \param inFileName name of input file (r-grid format)
      * \param outFileName name of output file (basis format)
      */
      void rGridToBasis(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Convert fields from Fourier (k-grid) to real-space (r-grid) format.
      *
      * \param inFileName name of input file (k-grid format)
      * \param outFileName name of output file (r-grid format)
      */
      void kGridToRGrid(const std::string& inFileName,
                        const std::string& outFileName);

      /**
      * Convert fields from real-space (r-grid) to Fourier (k-grid) format.
      *
      * \param inFileName name of input file (r-grid format)
      * \param outFileName name of output file (k-grid format)
      */
      void rGridToKGrid(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Convert fields from Fourier (k-grid) to symmetrized basis format.
      *
      * \param inFileName name of input file (k-grid format)
      * \param outFileName name of output file (basis format)
      */
      void kGridToBasis(const std::string& inFileName,
                        const std::string& outFileName);

      /**
      * Convert fields from symmetrized basis to Fourier (k-grid) format.
      *
      * \param inFileName name of input file (basis format)
      * \param outFileName name of output file (k-grid format)
      */
      void basisToKGrid(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Construct trial w-fields from c-fields.
      *
      * This function reads concentration fields in symmetrized basis
      * format and constructs an initial guess for corresponding chemical
      * potential fields by setting the Lagrange multiplier field xi to
      * zero. The resulting guess is stored in the System w fields 
      * container and is also output to a file in basis format.
      *
      * Upon return, w().hasData() and w().isSymmetric()are set true and
      * hasCFields is set false.
      *
      * \param inFileName  name of input c-field file (in, basis format)
      * \param outFileName  name of output w-field file (out, basis format)
      */
      void guessWfromC(const std::string& inFileName,
                       const std::string& outFileName);

      /**
      * Compare two field files in symmetrized basis format.
      *
      * Outputs maximum and root-mean-squared differences.
      *
      * \param field1  first array of fields (basis format)
      * \param field2  second array of fields (basis format)
      */
      void compare(const DArray< DArray<double> > field1,
                   const DArray< DArray<double> > field2);

      /**
      * Compare two field files in symmetrized basis format.
      *
      * Outputs maximum and root-mean-squared differences.
      *
      * \param field1  first array of fields (r-grid format)
      * \param field2  second array of fields (r-grid format)
      */
      void compare(const DArray< RField<D> > field1,
                   const DArray< RField<D> > field2);

      /**
      * Check if r-grid fields have the declared space group symmetry.
      *
      * \param inFileName name of input file
      * \return true if fields all have symmetry, false otherwise
      */
      bool checkRGridFieldSymmetry(const std::string & inFileName);

      //@}
      /// \name Member Accessors
      //@{

      /**
      * Get all of the chemical potential fields (const reference).
      */
      WFieldContainer<D> const & w() const;

      /**
      * Get all of the monomer concentration fields (const reference).
      */
      CFieldContainer<D> const & c() const;

      /**
      * Get all of the external potential fields (reference).
      */
      WFieldContainer<D>& h();

      /**
      * Get the mask (the field to which the total density is constrained) 
      * (reference).
      */
      Mask<D>& mask();

      /**
      * Get the Mixture by reference.
      */
      Mixture<D>& mixture();

      /**
      * Get Interaction (excess free energy model) by reference.
      */
      Interaction& interaction();

      /**
      * Get Interaction (excess free energy model) by const reference.
      */
      Interaction const & interaction() const;

      /**
      * Get Domain by const reference.
      */
      Domain<D> const & domain() const;

      /**
      * Get UnitCell (i.e., type and parameters) by const reference.
      */
      UnitCell<D> const & unitCell() const;

      /**
      * Get the spatial discretization mesh by const reference.
      */
      Mesh<D> const & mesh() const;

      /**
      * Get associated Basis object by reference.
      */
      Basis<D> const & basis() const;

      /**
      * Get associated FFT object.
      */
      FFT<D> const & fft() const;

      /**
      * Get associated FieldIo object.
      */
      FieldIo<D> const & fieldIo() const;

      /**
      * Get the Mixture by const reference.
      */
      Mixture<D> const & mixture() const;

      /**
      * Get the iterator by reference.
      */
      Iterator<D>& iterator();

      /**
      * Get the iterator by const reference.
      */
      Iterator<D> const & iterator() const;

      /**
      * Get the compressor by reference.
      */
      Compressor<D>& compressor();

      /**
      * Get homogeneous mixture (for reference calculations).
      */
      Homogeneous::Mixture& homogeneous();

      /**
      * Get const homogeneous mixture (for reference calculations).
      */
      Homogeneous::Mixture const & homogeneous() const;

      /**
      * Get McSimulator container for Monte Carlo moves.
      */
      McSimulator<D>& mcSimulator();

      /**
      * Get FileMaster by reference.
      */
      FileMaster& fileMaster();

      /**
      * Get FileMaster by const reference.
      */
      FileMaster const & fileMaster() const;

      /**
      * Get Random object by reference.
      */
      Random& random();

      /**
      * Get the group name string.
      */
      std::string groupName() const;

      /**
      * Does this system have an Iterator object?
      */
      bool hasIterator() const;

      /**
      * Does this system have a Sweep object?
      */
      bool hasSweep() const;

      /**
      * Does this system have external potential fields?
      */
      bool hasExternalFields() const;

      /**
       * Does this system have a mask (field to which the total density 
       * is constrained)?
       */
      bool hasMask() const;

      /**
      * Does this system have a Compressor object?
      */
      bool hasCompressor() const;

      /**
      * Does this system have an initialized McSimulator?
      */
      bool hasMcSimulator() const;

      /**
      * Have c fields been computed from the current w fields?
      */
      bool hasCFields() const;

      /**
      * Has the free energy been computed from the current w fields?
      */
      bool hasFreeEnergy() const;

      ///@}

   private:

      // Private member variables

      /**
      * Mixture object (solves MDE for all species).
      */
      Mixture<D> mixture_;

      /**
      * Domain object (unit cell, space group, mesh, and basis).
      */
      Domain<D> domain_;

      /**
      * Container for McMove objects (Monte Carlo Moves).
      */
      McSimulator<D> mcSimulator_;

      /**
      * Filemaster (holds paths to associated I/O files).
      */
      FileMaster fileMaster_;

      /**
      * Random number generator.
      */
      Random random_;

      /**
      * Homogeneous mixture, for reference.
      */
      Homogeneous::Mixture homogeneous_;

      /**
      * Pointer to Interaction (free energy model).
      */
      Interaction* interactionPtr_;

      /**
      * Pointer to an iterator.
      */
      Iterator<D>* iteratorPtr_;

      /**
      * Pointer to iterator factory object
      */
      IteratorFactory<D>* iteratorFactoryPtr_;

      /**
      * Pointer to a Sweep object
      */
      Sweep<D>* sweepPtr_;

      /**
      * Pointer to SweepFactory object
      */
      SweepFactory<D>* sweepFactoryPtr_;

      /**
      * Pointer to an compressor.
      */
      Compressor<D>* compressorPtr_;

      /**
      * Pointer to compressor factory object
      */
      CompressorFactory<D>* compressorFactoryPtr_;

      /**
      * Chemical potential fields.
      */
      WFieldContainer<D> w_;

      /**
      * Monomer concentration / volume fraction fields.
      */
      CFieldContainer<D> c_;

      /**
      * External potential fields.
      */
      WFieldContainer<D> h_;

      /**
      * Field to which the total density is constrained.
      */
      Mask<D> mask_;

      /**
      * Work array of field coefficients for all monomer types.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      mutable DArray< DArray<double> > tmpFieldsBasis_;

      /**
      * Work array of fields on real space grid.
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      mutable DArray< RField<D> > tmpFieldsRGrid_;

      /**
      * Work array of fields on Fourier grid (k-grid).
      *
      * Indexed by monomer typeId, size = nMonomer.
      */
      mutable DArray<RFieldDft<D> > tmpFieldsKGrid_;

      /**
      * Helmholtz free energy per monomer / kT.
      */
      double fHelmholtz_;

      /**
      * Ideal gas contribution to fHelmholtz_. 
      * 
      * This encompasses the internal energy and entropy of 
      * non-interacting free chains in their corresponding 
      * potential fields defined by w_.
      */
      double fIdeal_;

      /**
      * Multi-chain interaction contribution to fHelmholtz_.
      */
      double fInter_;

      /**
      * External field contribution to fHelmholtz_.
      */
      double fExt_;

      /**
      * Pressure times monomer volume / kT.
      */
      double pressure_;

      /**
      * Has the mixture been initialized?
      */
      bool hasMixture_;

      /**
      * Has the McSimulator been initialized?
      */
      bool hasMcSimulator_;

      /**
      * Has memory been allocated for fields in grid format?
      */
      bool isAllocatedRGrid_;

      /**
      * Has memory been allocated for fields in symmetrized basis format?
      */
      bool isAllocatedBasis_;

      /**
      * Have c fields been computed for the current w fields?
      *
      * Set true when c fields are computed by solving the MDEs for
      * all blocks, and set false whenever w fields or the unit cell
      * parameters are reset. When hasCFields_ is true, both the 
      * c fields for individual blocks and solvent species in the
      * Mixture and the fields for different monomer types the 
      * System::c_ container are those obtained from the current w 
      * fields in System::w_ container.
      */
      bool hasCFields_;

      /**
      * Has the free energy been computed for the current w and c fields?
      */ 
      bool hasFreeEnergy_;

      /**
      * Has the MC Hamiltonian been computed for the current w and c fields?
      */ 
      bool hasMcHamiltonian_;

      // Private member functions

      /**
      * Allocate memory for fields in grid formats (private)
      */
      void allocateFieldsGrid();

      /**
      * Allocate memory for fields in basis format (private)
      */
      void allocateFieldsBasis();

      /**
      * Read field file header and allocate memory for basis fields.
      *
      * \param filename name of field file
      */
      void allocateFieldsBasis(std::string filename);

      /**
      * Initialize Homogeneous::Mixture object.
      */
      void initHomogeneous();

      /**
      * Read a filename string and echo to log file.
      *
      * Used to read filenames in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param string  string to read and echo
      */
      void readEcho(std::istream& in, std::string& string) const;

   };

   // Inline member functions

   // Get the Domain<D> object.
   template <int D>
   inline Domain<D> const & System<D>::domain() const
   { return domain_; }

   // Get the UnitCell<D> object.
   template <int D>
   inline UnitCell<D> const & System<D>::unitCell() const
   { return domain_.unitCell(); }

   // Get the Mesh<D> object.
   template <int D>
   inline Mesh<D> const & System<D>::mesh() const
   { return domain_.mesh(); }

   // Get the associated Mixture object.
   template <int D>
   inline Mixture<D>& System<D>::mixture()
   { return mixture_; }

   // Get the associated const Mixture object.
   template <int D>
   inline Mixture<D> const & System<D>::mixture() const
   { return mixture_; }

   // Get the Basis<D> object.
   template <int D>
   inline Basis<D> const & System<D>::basis() const
   {  return domain_.basis(); }

   // Get the FFT<D> object.
   template <int D>
   inline FFT<D> const & System<D>::fft() const
   {  return domain_.fft(); }

   // Get the FieldIo<D> object.
   template <int D>
   inline FieldIo<D> const & System<D>::fieldIo() const
   {  return domain_.fieldIo(); }

   // Get the groupName string.
   template <int D>
   inline std::string System<D>::groupName() const
   { return domain_.groupName(); }

   // Get the McSimulator.
   template <int D>
   inline McSimulator<D>& System<D>::mcSimulator()
   {  return mcSimulator_; }

   // Get the FileMaster.
   template <int D>
   inline FileMaster& System<D>::fileMaster()
   {  return fileMaster_; }

   // Get the const FileMaster.
   template <int D>
   inline FileMaster const & System<D>::fileMaster() const
   {  return fileMaster_; }

   // Get the Random object.
   template <int D>
   inline Random& System<D>::random()
   {  return random_; }

   // Get the Homogeneous::Mixture object.
   template <int D>
   inline Homogeneous::Mixture& System<D>::homogeneous()
   {  return homogeneous_; }

   // Get the const Homogeneous::Mixture object.
   template <int D>
   inline Homogeneous::Mixture const & System<D>::homogeneous() const
   {  return homogeneous_; }

   // Get the Interaction (excess free energy model).
   template <int D>
   inline Interaction& System<D>::interaction()
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the const Interaction (excess free energy model).
   template <int D>
   inline Interaction const & System<D>::interaction() const
   {
      UTIL_ASSERT(interactionPtr_);
      return *interactionPtr_;
   }

   // Get the Iterator.
   template <int D>
   inline Iterator<D>& System<D>::iterator()
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   // Get the const Iterator.
   template <int D>
   inline Iterator<D> const & System<D>::iterator() const
   {
      UTIL_ASSERT(iteratorPtr_);
      return *iteratorPtr_;
   }

   // Get the Compressor
   template <int D>
   inline Compressor<D>& System<D>::compressor()
   {
      UTIL_ASSERT(compressorPtr_);
      return *compressorPtr_;
   }

   // Get container of chemical potential fields (const reference)
   template <int D>
   inline
   WFieldContainer<D> const & System<D>::w() const
   {  return w_; }

   // Get container of monomer concentration fields (const reference)
   template <int D>
   inline
   CFieldContainer<D> const & System<D>::c() const
   {  return c_; }

   // Get container of external potential fields (reference)
   template <int D>
   inline WFieldContainer<D>& System<D>::h()
   {  return h_; }

   // Get mask field (reference)
   template <int D>
   inline Mask<D>& System<D>::mask()
   {  return mask_; }

   // Does the system have an Iterator object?
   template <int D>
   inline bool System<D>::hasIterator() const
   {  return (iteratorPtr_ != 0); }

   // Does the system have a Sweep object?
   template <int D>
   inline bool System<D>::hasSweep() const
   {  return (sweepPtr_ != 0); }

   // Does this system have external potential fields? 
   template <int D>
   inline bool System<D>::hasExternalFields() const
   {  return h_.hasData(); }

   // Does this system have a mask?
   template <int D>
   inline bool System<D>::hasMask() const
   {  return mask_.hasData(); }

   // Does the system have a Compressor object?
   template <int D>
   inline bool System<D>::hasCompressor() const
   {  return (compressorPtr_ != 0); }

   // Does the system have an initialized McSimulator ?
   template <int D>
   inline bool System<D>::hasMcSimulator() const
   {  return (hasMcSimulator_); }

   // Have the c fields been computed for the current w fields?
   template <int D>
   inline bool System<D>::hasCFields() const
   {  return hasCFields_; }

   // Get the precomputed Helmholtz free energy per monomer / kT.
   template <int D>
   inline double System<D>::fHelmholtz() const
   {
      UTIL_CHECK(hasFreeEnergy_);
      return fHelmholtz_; 
   }

   // Get the precomputed pressure (units of kT / monomer volume).
   template <int D>
   inline double System<D>::pressure() const
   {  
      UTIL_CHECK(hasFreeEnergy_);
      return pressure_; 
   }

   // Has the free energy been computed for the current w fields?
   template <int D>
   inline bool System<D>::hasFreeEnergy() const
   {  return hasFreeEnergy_; }

   #ifndef PSPC_SYSTEM_TPP
   // Suppress implicit instantiation
   extern template class System<1>;
   extern template class System<2>;
   extern template class System<3>;
   #endif

} // namespace Pspc
} // namespace Pscf
#endif
