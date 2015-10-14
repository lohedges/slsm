/*
 * File:	InputOutput.h
 * Author:	lh627
 */

#ifndef _INPUTOUTPUT_H
#define _INPUTOUTPUT_H

class InputOutput
{
public:
	InputOutput();

    //! Save the level set function as a ParaView VTK file.
    /*! 
     */
    void saveLevelSetVTK();
};

#endif	/* _INPUTOUTPUT_H */
