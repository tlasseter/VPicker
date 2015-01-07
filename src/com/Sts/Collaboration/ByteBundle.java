/*
 * Created on Aug 11, 2004
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Generation - Code and Comments
 */
package com.Sts.Collaboration;

import java.io.*;


/**
 * 
 * 
 * @author druths
 */
public class ByteBundle implements Serializable {

  public boolean hasInt = false;
  public int val;
  public byte[] array;
  
  // constructors
  public ByteBundle(int val) {
    this.val = val;
    hasInt = true;
  }
  
  public ByteBundle(byte[] array) {
    this.array = array;
  }
}
