#ifndef RPC_SYSTEM_H
#define RPC_SYSTEM_H

/*
* PSCF - Polymer Self-Consistent Field Theory
*
* Copyright 2016 - 2022, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>     // base class

#include <rpc/solvers/Mixture.h>          // member
#include <rpc/field/Domain.h>             // member
#include <rpc/field/FieldIo.h>            // member
#include <rpc/field/WFieldContainer.h>    // member
#include <rpc/field/CFieldContainer.h>    // member
#include <rpc/field/Mask.h>               // member

#include <prdc/cpu/RField.h>               // member
#include <prdc/cpu/RFieldDft.h>            // member

#include <pscf/homogeneous/Mixture.h>      // member

#include <util/misc/FileMaster.h>          // member
#include <util/containers/DArray.h>        // member template
#include <util/containers/FSArray.h>       // member template

namespace Pscf {
   class Interaction;
   namespace Prdc {
      template <int D> class UnitCell;
   }
}

namespace Pscf {
namespace Rpc {

   // Forward references in Pscf::Rpc
   template <int D> class Iterator;
   template <int D> class IteratorFactory;
   template <int D> class Sweep;
   template <int D> class SweepFactory;
   template <int D> class Compressor;
   template <int D> class CompressorFactory;
   template <int D> class Simulator;
   template <int D> class SimulatorFactory;

   using namespace Util;
   using namespace Prdc;
   using namespace Prdc::Cpu;

   /**
   * Main class for SCFT simulation of one system.
   *
   * A System has (among other components):
   *
   *    - a Mixture (a container for polymer and solvent solvers)
   *    - an Interaction (list of binary chi parameters)
   *    - a Domain (description of unit cell and discretization)
   *    - a container of monomer chemical potential fields 
   *    - a container of monomer concentration fields 
   *
   * A system may also optionally contain Iterator, Sweep, Compressor
   * and Simulator objects.
   *
   * Typical usage of a System<D> object looks something like this:
   * \code
   *    System<D> system;
   *    system.setOptions(argc, argv);
   *    system.readParam();
   *    system.readCommands();
   * \endcode
   * where argc, and argv are parameters containing information about
   * command line arguments that must be passed from the main program.
   * This is implemented as function template Pscf::Rpc::run in the
   * file src/rpc/pscf_pc.cpp.
   *
   * \ref user_param_pc_page "Parameter File Format"
   * \ref user_command_pc_page "Command File Format"
   * \ingroup Pscf_Rpc_Module
   */
   template <int D>
   class System : public ParamComposite
   {

   public:

      /// \name Construction and Destruction
      ///@{

      /**
      * Constructor.
      */
      System();

      /**
      * Destructor.
      */
      ~System();

      ///@}
      /// \name Lifetime (Actions)
      ///@{
         
      double computeDebye(double x);
      double computeIntraCorrelation(double qSquare);
      RField<D> computeIntraCorrelations();
      RField<D> estimateLR(double stepSize);
      void estimateError(double stepSize);
      double a_;
      double f_;

      /**
      * Process command line options.
      *
      * This function takes the same arguments as any C/C++ main program
      * function. The arguments of the main function should d be passed 
      * to this function unaltered, to allow this function to process the
      * command line options.
      *
      * \param argc number of command line arguments
      * \param argv array of pointers to command line arguments
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
      *
      * This function reads the parameter file set by the -p command
      * line option.
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
      *
      * This function reads the parameter file set by the -c command
      * line option.
      */
      void readCommands();

      ///@}
      /// \name W Field Modifiers
      ///@{

      /**
      * Read chemical potential fields in symmetry adapted basis format.
      *
      * This function opens and reads the file with the name given by the
      * "filename" string parameter, which must contain chemical potential 
      * fields in symmetry-adapted basis format. The function sets the
      * system w fields equal to those given in this file, by copying
      * elements of the representation in basis format and computing the
      * representation in r-grid format. On exit, both w().basis() and
      * w().rgrid() have been reset, w().hasData and w().isSymmetric() 
      * are true, and hasCFields() is false.
      *
      * \param filename name of input w-field basis file
      */
      void readWBasis(const std::string & filename);

      /**
      * Read chemical potential fields in real space grid (r-grid) format.
      *
      * This function opens and reads the file with the name given by the
      * "filename" string, which must contain chemical potential fields
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
      * not set components the symmetry-adapted basis format. On return,
      * w.rgrid() is reset, w().hasData() is true, w().isSymmetric() is 
      * false, and hasCFields() is false.
      *
      * \param fields  array of new w (chemical potential) fields
      */
      void setWRGrid(DArray< RField<D> > const & fields);

      /**
      * Construct trial w-fields from c-fields.
      *
      * This function reads concentration fields in symmetrized basis
      * format and constructs an initial guess for corresponding chemical
      * potential fields by setting the Lagrange multiplier field xi to
      * zero. The result is stored in the System w fields container
      *
      * Upon return, w().hasData() and w().isSymmetric() are set true,
      * while hasCFields is set false.
      *
      * \param filename  name of input c-field file (basis format)
      */
      void estimateWfromC(const std::string& filename);

      ///@}
      /// \name Unit Cell Modifiers
      ///@{

      /**
      * Set parameters of the associated unit cell.
      *
      * The lattice set in this UnitCell must agree with any lattice
      * value that was set previously in the parameter file.
      * 
      * \param unitCell  new UnitCell<D> (i.e., new parameters)
      */
      void setUnitCell(UnitCell<D> const & unitCell);

      /**
      * Set state of the associated unit cell.
      *
      * The lattice parameter must agree with any lattice value that
      * was set previously in the parameter file.
      * 
      * \param lattice  lattice system
      * \param parameters  array of new unit cell parameters.
      */
      void setUnitCell(typename UnitCell<D>::LatticeSystem lattice,
                       FSArray<double, 6> const & parameters);

      /**
      * Set parameters of the associated unit cell.
      *
      * The size of the FSArray<double> parameters must match the
      * expected number of parameters for the current lattice type.
      * 
      * \param parameters  array of new unit cell parameters.
      */
      void setUnitCell(FSArray<double, 6> const & parameters);

      ///@}
      /// \name Primary Field Theory Computations
      ///@{

      /**
      * Solve the modified diffusion equation once, without iteration.
      *
      * This function calls the Mixture::compute() function to solve
      * the statistical mechanics problem for a non-interacting system
      * subjected to the currrent chemical potential fields. This
      * requires solution of the modified diffusion equation for all 
      * polymers, computation of Boltzmann weights for all solvents, 
      * computation of molecular partition functions for all species, 
      * computation of concentration fields for blocks and solvents, 
      * and computation of overall concentrations for all monomer types. 
      * This function does not compute the canonical (Helmholtz) free 
      * energy or grand-canonical free energy (i.e., pressure). Upon 
      * return, the flag hasCFields is set true.
      *
      * If argument needStress == true, then this function also calls
      * Mixture<D>::computeStress() to compute the stress.
      *
      * \pre The w().hasData() flag must be true on entry, to confirm
      * that chemical potential fields have been set.
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
      * \pre The w().hasData() flag must be true on entry, to confirm
      * that chemical potential fields have been set. 
      *
      * \pre The w().isSymmetric() flag must be set true if the chosen 
      * iterator uses a basis representation, and thus requires this.
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
      
      /**
      * Write timer file to an ostream
      *
      * \param out output stream
      */
      void writeTimers(std::ostream& out);
      
      /**
      * Clear timers 
      */
      void clearTimers();
      
      ///@}
      /// \name Thermodynamic Properties
      ///@{

      /**
      * Compute free energy density and pressure for current fields.
      *
      * This function should be called after a successful call of
      * System::iterate() or Iterator::solve(). Resulting values are 
      * stored and then accessed by the fHelmholtz() and pressure() 
      * functions.
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

      ///@}
      /// \name Thermodynamic Data Output 
      ///@{

      /**
      * Write parameter file to an ostream, omitting any sweep block.
      *
      * This function omits the Sweep block of the parameter file, if
      * any, in order to allow the output produced during a sweep to refer
      * only to parameters relevant to a single state point, and to be
      * rerunnable as a parameter file for a single SCFT calculation.
      *
      * \param out output stream
      */
      void writeParamNoSweep(std::ostream& out) const;

      /**
      * Write thermodynamic properties to a file.
      *
      * This function outputs Helmholtz free energy per monomer, pressure
      * (in units of kT per monomer volume), the volume fraction and 
      * chemical potential of each species, and all unit cell parameters.
      * 
      * If parameter "out" is a file that already exists, this function
      * will append this information to the end of the file, rather than 
      * overwriting that file. Calling writeParamNoSweep and writeThermo 
      * in succession with the same file will thus produce a single file
      * containing both input parameters and resulting thermodynanic
      * properties.
      *
      * \param out output stream
      */
      void writeThermo(std::ostream& out);

      ///@}
      /// \name Field Output 
      ///@{

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
      * appear ordered by polymer id and then by block id, with blocks
      * of the same polymer listed sequentially, followed by columns
      * associated with solvent species ordered by solvent id.
      *
      * \param filename name of output file
      */
      void writeBlockCRGrid(const std::string & filename) const;

      ///@}
      /// \name Propagator Output 
      ///@{

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

      ///@}
      /// \name Crystallographic Information
      ///@{

      /**
      * Output information about stars and symmetrized basis functions.
      *
      * This function opens a file with the specified filename, calls
      * Basis<D>::outputStars, and closes the file before returning.
      *
      * \param filename name of output file 
      */
      void writeStars(std::string const & filename) const;

      /**
      * Output information about waves.
      *
      * This function opens a file with the specified filename, calls
      * Basis<D>::outputWaves, and closes the file before returning.
      *
      * \param filename name of output file 
      */
      void writeWaves(std::string const & filename) const;

      /**
      * Output all elements of the space group.
      *
      * \param filename name of output file 
      */
      void writeGroup(std::string const & filename) const;

      ///@}
      /// \name Field File Manipulations 
      ///@{

      /**
      * Convert a field from symmetrized basis format to r-grid format.
      *
      * This function reads a field file in basis format, converts the 
      * fields to r-grid format, and writes the fields in r-grid format 
      * to a different file. 
      *
      * This and other field conversion functions do not change the w 
      * or c fields stored by this System - all required calculations 
      * are performed using temporary or mutable memory. 
      * 
      * \param inFileName name of input file (basis format)
      * \param outFileName name of output file (r-grid format)
      */
      void basisToRGrid(const std::string & inFileName,
                        const std::string & outFileName);

      /**
      * Convert a field from real-space grid to symmetrized basis format.
      *
      * This function checks if the input fields have the declared space
      * group symmetry, and prints a warning if it detects deviations
      * that exceed some small threshhold, but proceeds to attempt the
      * conversion even if such an error is detected. Converting a field 
      * that does not have the declared space group symmetry to basis
      * format is a destructive operation that modifies the field in
      * unpredictable ways.
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
      * This function checks if the input fields have the declared space
      * group symmetry, and prints a warning if it detects deviations
      * that exceed some small threshhold, but proceeds to attempt the
      * conversion even if such an error is detected. Converting a field 
      * that does not have the declared space group symmetry to basis
      * format is a destructive operation that modifies the field in
      * unpredictable ways.
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
      * \param epsilon error threshold used when testing for symmetry
      * \return true if fields all have symmetry, false otherwise
      */
      bool checkRGridFieldSymmetry(const std::string & inFileName,
                                   double epsilon = 1.0E-8);

      /**
      * Expand the number of spatial dimensions of an r-grid field.
      *
      * This function reads a D-dimensional field and outputs a field
      * in a format appropriate for d-dimensional space, for d > D, by
      * assuming that all field values are independent of coordinates
      * associated with the added dimensions. It can thus create a file
      * representing a field with lamellar (D=1) or hexagonal (D=2)
      * symmetry on a 3D (d=3) grid. 
      *
      * Element i of array newGridDimensions contains the number of
      * grid points in added dimension D + i. This array must have a
      * capacity d - D. 
      *
      * \param inFileName filename name of input field file
      * \param outFileName filename name of output field file
      * \param d  intended dimensions (d > D)
      * \param newGridDimensions number of grid points in added dimensions
      */
      void expandRGridDimension(const std::string & inFileName,
                                const std::string & outFileName,
                                int d,
                                DArray<int> newGridDimensions);
      /**
      * Replicate unit cell
      *
      * This function reads a D-dimensional field and replicate 
      * the unit cell a specified number of times in each D direction
      * 
      * Element i of array replicas contains the number of replication 
      * times in direction i. 
      * 
      * \param inFileName filename name of input field file
      * \param outFileName filename name of output field file
      * \param replicas  the number of replicas in each D direction
      */
      void replicateUnitCell(const std::string & inFileName,
                             const std::string & outFileName,
                             IntVec<D> const & replicas);

      ///@}
      /// \name Field Accessors
      ///@{

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
      * Get the mask (field to which total density is constrained). 
      */
      Mask<D>& mask();

      ///@}
      /// \name Member Object Accessors
      ///@{

      /**
      * Get the Mixture by non-const reference.
      */
      Mixture<D>& mixture();

      /**
      * Get the Mixture by const reference.
      */
      Mixture<D> const & mixture() const;

      /**
      * Get Interaction (excess free energy model) by reference.
      */
      Interaction& interaction();

      /**
      * Get Interaction (excess free energy model) by const reference.
      */
      Interaction const & interaction() const;

      /**
      * Get Domain by non-const reference.
      */
      Domain<D> const & domain() const;
      
      /**
      * Get Domain by const reference.
      */
      Domain<D> & domain();

      /**
      * Get UnitCell (i.e., type and parameters) by const reference.
      */
      UnitCell<D> const & unitCell() const;
      
      /**
      * Get UnitCell (i.e., type and parameters) by non-const reference.
      */
      UnitCell<D> & unitCell();

      /**
      * Get the spatial discretization mesh by const reference.
      */
      Mesh<D> const & mesh() const;

      /**
      * Get the Basis by const reference.
      */
      Basis<D> const & basis() const;

      /**
      * Get associated FFT object by const reference.
      */
      FFT<D> const & fft() const;

      /**
      * Get associated FieldIo object by const reference.
      */
      FieldIo<D> const & fieldIo() const;

      /**
      * Get the iterator.
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
      * Get Simulator for field theoretic simulation. 
      */
      Simulator<D>& simulator();

      /**
      * Get homogeneous mixture (for reference calculations).
      */
      Homogeneous::Mixture& homogeneous();

      /**
      * Get const homogeneous mixture (for reference calculations).
      */
      Homogeneous::Mixture const & homogeneous() const;

      /**
      * Get the FileMaster.
      *
      * Access by non-const reference is used in some unit tests.
      */
      FileMaster& fileMaster();

      /**
      * Get the FileMaster by const reference.
      */
      FileMaster const & fileMaster() const;

      /**
      * Get the group name string.
      */
      std::string groupName() const;

      ///@}
      /// \name Queries
      ///@{

      /**
      * Does this system have an Iterator object?
      */
      bool hasIterator() const;

      /**
      * Have c fields been computed from the current w fields?
      */
      bool hasCFields() const;

      /**
      * Has the free energy been computed from the current w fields?
      */
      bool hasFreeEnergy() const;

      /**
      * Does this system have external potential fields?
      */
      bool hasExternalFields() const;

      /**
      * Does this system have a mask (inhomogeneous density constraint)
      */
      bool hasMask() const;

      /**
      * Does this system have a Sweep object?
      */
      bool hasSweep() const;

      /**
      * Does this system have a Compressor object?
      */
      bool hasCompressor() const;

      /**
      * Does this system have an initialized Simulator?
      */
      bool hasSimulator() const;

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
      * Filemaster (holds paths to associated I/O files).
      */
      FileMaster fileMaster_;

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
      * Simulator - coordinator for field theoretic simulation.
      */
      Simulator<D>* simulatorPtr_;

      /**
      * Pointer to Simulator factory object
      */
      SimulatorFactory<D>* simulatorFactoryPtr_;

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
      mutable DArray< RFieldDft<D> > tmpFieldsKGrid_;

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
      * 
      * This is -1 times the grand-canonical free energy per monomer,
      * divided by kT.
      */
      double pressure_;

      /**
      * Has the mixture been initialized?
      */
      bool hasMixture_;

      /**
      * Has memory been allocated for fields in grid format?
      */
      bool isAllocatedGrid_;

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
      * Has fHelmholtz been computed for the current w and c fields?
      */ 
      bool hasFreeEnergy_;

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
      * Read a field file header, make the basis if not done previously.
      * 
      * Used to peek at a file header to get initial unit cell parameters,
      * use this to initialize basis if not done previously.
      *
      * \param filename name of field file
      */
      void readFieldHeader(std::string filename);

      /**
      * Read a string and echo to log file.
      *
      * Used to read filenames in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param string  string to read and echo
      */
      void readEcho(std::istream& in, std::string& string) const;

      /**
      * Read a floating point number and echo to log file.
      *
      * Used to read filenames in readCommands.
      *
      * \param in  input stream (i.e., input file)
      * \param value  number to read and echo
      */
      void readEcho(std::istream& in, double& value) const;

      /**
      * Initialize Homogeneous::Mixture object.
      */
      void initHomogeneous();

   };

   // Inline member functions

   // Get the Mixture object.
   template <int D>
   inline Mixture<D>& System<D>::mixture()
   { return mixture_; }

   // Get the Mixture by const reference.
   template <int D>
   inline Mixture<D> const & System<D>::mixture() const
   { return mixture_; }

   // Get the Domain by const reference.
   template <int D>
   inline Domain<D> const & System<D>::domain() const
   { return domain_; }
   
   // Get the Domain by nonconst reference.
   template <int D>
   inline Domain<D> & System<D>::domain() 
   { return domain_; }

   // Get the UnitCell by const reference.
   template <int D>
   inline UnitCell<D> const & System<D>::unitCell() const
   { return domain_.unitCell(); }
   
   // Get the UnitCell by const reference.
   template <int D>
   inline UnitCell<D> & System<D>::unitCell()
   { return domain_.unitCell(); }

   // Get the Mesh by const reference.
   template <int D>
   inline Mesh<D> const & System<D>::mesh() const
   { return domain_.mesh(); }

   // Get the Basis<D> object.
   template <int D>
   inline Basis<D> const & System<D>::basis() const
   {  return domain_.basis(); }

   // Get the FFT object by const reference.
   template <int D>
   inline FFT<D> const & System<D>::fft() const
   {  return domain_.fft(); }

   // Get the FieldIo<D> by const reference.
   template <int D>
   inline FieldIo<D> const & System<D>::fieldIo() const
   {  return domain_.fieldIo(); }

   // Get the groupName string by value.
   template <int D>
   inline std::string System<D>::groupName() const
   { return domain_.groupName(); }

   // Get the Simulator.
   template <int D>
   inline Simulator<D>& System<D>::simulator()
   {  return *simulatorPtr_; }

   // Get the FileMaster by non-const reference.
   template <int D>
   inline FileMaster& System<D>::fileMaster()
   {  return fileMaster_; }

   // Get the FileMaster by const reference.
   template <int D>
   inline FileMaster const & System<D>::fileMaster() const
   {  return fileMaster_; }

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

   // Get the Interaction by const reference.
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

   // Get the Iterator by const reference.
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

   // Does the system have an initialized Simulator ?
   template <int D>
   inline bool System<D>::hasSimulator() const
   {  return (simulatorPtr_ != 0); }

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

   #ifndef RPC_SYSTEM_TPP
   // Suppress implicit instantiation
   extern template class System<1>;
   extern template class System<2>;
   extern template class System<3>;
   #endif

} // namespace Rpc
} // namespace Pscf
#endif
