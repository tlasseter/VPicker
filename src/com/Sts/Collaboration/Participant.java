package com.Sts.Collaboration;

import java.io.*;
import java.net.*;

/**
 * This class contains all information pertinent to a given participant.
 * 
 * @author druths
 */
public class Participant implements Serializable {
	// fields
	private InetAddress _addr;
	private int _port;

    // constructors
	public Participant(InetAddress addr, int port) {
		_addr = addr;
		_port = port;
	}
	
    //// Methods
	public InetAddress getAddress() {
		return _addr;
	}
	
    public int getPort() {
		return _port;
	}
    
    public String toString() {
    	    return "" + _addr + ":" + _port;
    }
    
    public boolean equals(Object obj) {
     
        if(!(obj instanceof Participant)) {
         
            return false;
        } else {
         
            Participant p2 = (Participant) obj;
            
            return (p2._addr.equals(_addr) && p2._port == _port);
        }
    }
    
    public int hashCode() {
    	    return _port;
    }
}