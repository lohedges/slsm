/*
 * File:	InputOutput.h
 * Author:	lester
 */

#ifndef _INPUTOUTPUT_H
#define _INPUTOUTPUT_H

#include <fstream>
#include <sstream>
#include <string>

#include "Boundary.h"

class InputOutput
{
public:
    //! Constructor.
	InputOutput();

    //! Save the level set function as a ParaView VTK file.
    /*! \param datapoint
            The datapoint of the current optimisation trajectory.

        \param mesh
            A reference to the finite element mesh.

        \param levelSet
            A reference to the level set object.

        \param outputDirectory
            The output directory path (optional).
     */
    void saveLevelSetVTK(const unsigned int&, const Mesh&, const LevelSet&, const std::string& outputDirectory = "") const;

    //! Save the level set function as a ParaView VTK file.
    /*! \param fileName
            The name of the data file.

        \param mesh
            A reference to the finite element mesh.

        \param levelSet
            A reference to the level set object.
     */
    void saveLevelSetVTK(const std::ostringstream&, const Mesh&, const LevelSet&) const;

    //! Save the level set function as a plain text file.
    /*! \param datapoint
            The datapoint of the current optimisation trajectory.

        \param mesh
            A reference to the finite element mesh.

        \param levelSet
            A reference to the level set object.

        \param outputDirectory
            The output directory path (optional).

        \param isXY
            Whether to also output the nodal x/y coordinates (optional).
     */
    void saveLevelSetTXT(const unsigned int&, const Mesh&, const LevelSet&,
        const std::string& outputDirectory = "", bool isXY = false) const;

    //! Save the level set function as a plain text file.
    /*! \param fileName
            The name of the data file.

        \param mesh
            A reference to the finite element mesh.

        \param levelSet
            A reference to the level set object.

        \param isXY
            Whether to also output the nodal x/y coordinates (optional).
     */
    void saveLevelSetTXT(const std::ostringstream&, const Mesh&, const LevelSet&, bool isXY = false) const;
};

#endif	/* _INPUTOUTPUT_H */
