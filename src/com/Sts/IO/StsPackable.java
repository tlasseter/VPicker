
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.IO;

import com.Sts.Utilities.*;

public interface StsPackable
{
  	public String packFields() throws StsException;
  	public void unpackFields(String s) throws StsException;
}
