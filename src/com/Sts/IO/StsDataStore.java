
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.IO;

import com.Sts.Utilities.*;

import java.io.*;

abstract public class StsDataStore
{
    // file access status
    final static public int NO_ACCESS = 0;
    final static public int WRITE_ACCESS = 1;
    final static public int READ_ACCESS = 2;

    private String name = null;  // data store name (e.g. filename)
    private int access = NO_ACCESS;

    /** default constructor */
    public StsDataStore() { }


    // abstract methods callable only from descendents

    /** open this data store for read or write */
    public abstract boolean openDataStore();

    /** close this data store */
    public abstract boolean closeDataStore();

    /** remove this data store */
    public abstract boolean removeDataStore(String name);

    // abstract methods callable from anywhere

    /** set/get a header string in this data store.
        this may be encoded metadata that must be parsed
        in the derived data store */
    public abstract void setHeader(String string);
    public abstract String getHeader();

    /** set/get a float vector in this data store */
//    public abstract boolean setFloatVector(float[] vector);
//    public abstract float[] getFloatVector();

    /** set/get a byte vector in this data store */
//    public abstract boolean setByteVector(byte[] vector);
//    public abstract byte[] getByteVector();


    // instance methods

    /** accessors */
    public String getName() { return name; }
    public boolean isOpen() { return name!=null; }
    public void setAccess(int access) { this.access = access; }
    public int getAccess() { return access; }
    public boolean readAccess() { return access==READ_ACCESS; }
    public boolean writeAccess() { return access==WRITE_ACCESS; }

    /** open this data store */
    public void open(String name, int access) throws StsException
    {
        if (this.name != null) closeDataStore();
        setAccess(access);
        this.name = name;
        if (name != null) openDataStore();
    }

    public void open(File rootDir, String name, int access) throws StsException
    {
        if (this.name != null) closeDataStore();
        setAccess(access);
        this.name = rootDir.getAbsolutePath() + File.separator + name;
        if (name != null) openDataStore();
    }

    /** close this data store */
    public void close() throws StsException
    {
        if (name == null) return;

        closeDataStore();
        name = null;
    }

    /** remove this data store */
    public void remove() throws StsException
    {
        if (name==null) return;

        closeDataStore();  // close it first
        removeDataStore(name);
        name = null;
    }

    /** remove any data store we have the name for */
    public void remove(String name) throws StsException
    {
        if (name==null) return;

        if (name.equals(this.name)) remove();  // this data store
        else removeDataStore(name);
    }


}
