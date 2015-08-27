//
//  individual.cpp
//  SHERIF
//
//  Created by David CHAMPREDON on 2015-02-25.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#include "individual.h"



void individual::create()
{
	_maxstate	= 99999; // <-- change that to 2(nE+nI)+nH+nF+3 = getState_D() ; but need to be safely coded...

	// susceptible in general population by default
	_state		= 0;
	_isHCW		= false;
	
	_timeDiseaseAcquisition	= 9E9;	// very large time means not infected yet
	_infectorID				= 0;
	_GIbck					= -999.999;
	_GIfwd.resize(0);
}



void individual::incrementStatus()
{
	if(_state<_maxstate) _state++;
}


void individual::displayInfo()
{
	cout<<endl;
	cout<<"Infectious status = " << _state << endl;
	cout<<"Infection time = " << _timeDiseaseAcquisition << endl;
	cout<<"Infector ID = " << _infectorID << endl;
	cout<<"GIbck = " << _GIbck<< endl;
	cout<<"GIfwd:";
	if(_GIfwd.size()>0) displayVector(_GIfwd);
	cout<<endl;
}