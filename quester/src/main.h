/**@page quester QUESTER: QUEST inTERpreter documentation
 *
 * @section questerintro Introduction
 *
 * QUESTER is a program that can run the QUEST simulation programs. It allows
 * both interactive (MATLAB like) and batch mode (python like) of operations. 
 *
 * Features:
 * - MATLAB like shell allows you to iteractively build and run simulations.
 *     - Allows visualization of structure and results.
 * - Python like batch mode allows you to run the simulations batch mode.
 *     - Allows users to run simulations in parallel in 
 *       cluster/supercomputers.
 *
 * @section questerexample Examples
 *
 * - Interactive mode:
 *
 * @code {SH}
 * $ quester
 *
 *                             Q  U  E  S  T  E  R
 *                                   vX.YY.ZZ
 *                                quest-vX.YY.ZZ

 * For help, type help
 *
 * |0> 
 *
 * @endcode
 *
 * - Interactive mode:
 *
 * @code {SH}
 * $ quester
 *
 *                             Q  U  E  S  T  E  R
 *                                   vX.YY.ZZ
 *                                quest-vX.YY.ZZ
 * For help, type help
 *
 * |0> graphene = CrystalStruct.createPrimitiveCell (CrystalStruct.GRAPHENE)
 *     <GrapheneCell>
 * |1> hamiltonian = Hamiltonian.createHamiltonian (graphene)
 *     <Hamiltonian>
 * |2> calculator = Calculator.BandStructure ()
 * |3> kpoints = BrillouinZone.getLines (graphene, BrillouinZone.G, 
 *         ... BrillouinZone.M, BrillouinZone.K, BrillouinZone.G, 100);
 * |4> Ek = calculator.run (kpoints)
 * |5> plotEk (Ek, kpoints)
 *
 * @endcode
 *
 * - Interactive mode -- run a script:
 *
 * @code {SH}
 * $ quester -i -r "hello.py"
 *
 *                             Q  U  E  S  T  E  R
 *                                   vX.YY.ZZ
 *                                quest-vX.YY.ZZ
 *
 * For help, type help
 *
 * |0> run hello.py
 *     Hello world!
 *
 * @endcode
 *
 * - Batch mode -- run a script:
 *
 * @code {SH}
 * $ quester -b -r "hello.py"
 *
 *                             Q  U  E  S  T  E  R
 *                                   vX.YY.ZZ
 *                                quest-vX.YY.ZZ
 *
 * Hello world!
 *
 * $
 *
 * @endcode

 *
 */

#ifndef QUESTER_MAIN_H
#define QUESTER_MAIN_H


#endif



