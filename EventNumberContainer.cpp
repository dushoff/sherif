//
//  EventNumberContainer.cpp
//  SHERIF
//
//  Created by David CHAMPREDON on 2015-08-07.
//  Copyright (c) 2015 David CHAMPREDON. All rights reserved.
//

#include "EventNumberContainer.h"



unsigned long EventNumberContainer::sumAll_infections(){

	unsigned long N = _event_label.size();
	unsigned long res = 0;
	
	for(int i=0; i<N; i++)
		if(is_in_string(_event_label[i], "infection")) res += _n_event[i];
	
	return(res);
}


void EventNumberContainer::display(){
	unsigned long N = _event_label.size();
	cout << "Event Label \t sub-Label \t n Events"<<endl;
	coutline(40);
	for(int i=0; i<N; i++){
		cout << _event_label[i]<<"\t"<<_event_sublabel[i]<<"\t"<<_n_event[i]<<endl;
	}
}


void EventNumberContainer::display_nonzeroEvents(){
	unsigned long N = _event_label.size();
	cout << "Event Label \t sub-Label \t n Events"<<endl;
	coutline(40);
	for(int i=0; i<N; i++){
		if(_n_event[i]>0) cout << _event_label[i]<<"\t"<<_event_sublabel[i]<<"\t"<<_n_event[i]<<endl;
	}
	coutline(40);
}

