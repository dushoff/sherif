//
//  individual.h
//  SHERIF
//
//  Created by David CHAMPREDON on 2015-02-25.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#ifndef __SHERIF__individual__
#define __SHERIF__individual__

#include <stdio.h>
#include "dcTools.h"

class individual
{
	unsigned long	_ID;
	
	unsigned int	_state;
	unsigned int	_maxstate;
	
	double			_timeDiseaseAcquisition;
	unsigned long	_infectorID;
	vector<double>	_timeDiseaseTransmit;
	
	double			_GIbck;
	vector<double>	_GIfwd;

	bool			_isHCW;
	
public:

	individual(){};
	
	void create();
	void clean();
	
	
	// ====== SET FUNCTIONS =======
	
	void set_ID(unsigned long i) {_ID=i;}
	
	void set_state(unsigned int s) {_state=s;}
	void set_maxstate(unsigned int s) {_maxstate=s;}
	void set_isHCW(bool x) {_isHCW=x;}
	
	void set_timeDiseaseAcquisition(double x) {_timeDiseaseAcquisition=x;}
	void set_timeDiseaseTransmit(double x) {_timeDiseaseTransmit.push_back(x);}
	
	void set_infectorID(unsigned long i) {_infectorID=i;}
	
	void set_GIbck(double x) {_GIbck = x;}
	void set_GIfwd_incr(double x) {_GIfwd.push_back(x);}
	
	
	// ====== GET FUNCTIONS =======

	unsigned int		get_state() {return _state;}
	unsigned long		get_ID() {return _ID;}
	unsigned long		get_infectorID() {return _infectorID;}
	double				get_timeDiseaseAcquisition(){return _timeDiseaseAcquisition;}
	vector<double>		get_timeDiseaseTransmit(){return _timeDiseaseTransmit;}
	bool				get_isHCW(){return _isHCW;}
	
	double				get_GIbck(){return _GIbck;}
	vector<double>		get_GIfwd(){return _GIfwd;}
	
	
	
	// ====== EPIDEMIC =======
	
	void incrementStatus();
	
	
	// ====== HELPERS =======
	
	void displayInfo();
	
};




#endif /* defined(__SHERIF__individual__) */
