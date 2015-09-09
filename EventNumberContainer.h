//
//  EventNumberContainer.h
//  SHERIF
//
//  Created by David CHAMPREDON on 2015-08-07.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#ifndef __SHERIF__EventNumberContainer__
#define __SHERIF__EventNumberContainer__

#include "dcTools.h"

using namespace std;

struct EventNumberContainer{
	
	/// Container for number of events per event when
	/// using a tau-leap algorithm
	
	vector<string> _event_label;
	vector<unsigned int> _event_sublabel;
	vector<unsigned long> _n_event;
	
	EventNumberContainer(){};
	EventNumberContainer(vector<string> event_label,
						 vector<unsigned int> event_sublabel,
						 vector<unsigned long> n_event){
		_event_label = event_label;
		_event_sublabel = event_sublabel;
		_n_event = n_event;
	};

	void display();
	void display_nonzeroEvents();
	unsigned long sumAll_infections();

};


#endif /* defined(__SHERIF__EventNumberContainer__) */
