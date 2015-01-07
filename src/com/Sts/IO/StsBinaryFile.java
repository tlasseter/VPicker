

//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.IO;

import com.Sts.Utilities.*;

import java.io.*;

public class StsBinaryFile
{
    // instance fields
    StsAbstractFile file = null;
    String filename = null;
    DataInputStream dis = null;
    DataOutputStream dos = null;
    boolean includeSizeInFile = true;

    public StsBinaryFile(StsAbstractFile file)
    {
        this.file = file;
        filename = file.getFilename();
    }

    public boolean openRead()
    {
        try
        {
            if(dis != null) close();
            InputStream is = file.getInputStream();
            BufferedInputStream bis = new BufferedInputStream(is);
            dis = new DataInputStream(bis);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsBinaryFile.openReadAndCheck() failed." +
                "Can't read: " + file.getFilename(), e, StsException.WARNING);
            return false;
        }
    }

    public boolean openWrite()
    {
        return openWrite(true);
    }

    public boolean openWrite(boolean append)
    {
        try
        {
            if(dos != null) close();
            OutputStream os = file.getOutputStream(append); // true: append write to end of file
            BufferedOutputStream bos = new BufferedOutputStream(os);
            dos = new DataOutputStream(bos);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsBinaryFile.openWrite() failed." +
                "Can't write: " + filename, e, StsException.WARNING);
            return false;
        }
    }

    public boolean openReadWrite()
    {
        return openRead() && openWrite();
    }

    /** close this binary file */
    public boolean close()
    {
        try
        {
            if (dos != null)
            {
                dos.flush();
                dos.close();
                dos = null;
            }
            if (dis != null)
            {
                dis.close();
                dis = null;
            }
			return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.close() failed."
                    + "Unable to close file " + filename, e, StsException.WARNING);
			return false;
        }
    }
    /** set/get an float values in this binary file */
    public boolean setFloatValues(float[] vector, boolean incSize)
    {
        includeSizeInFile = incSize;
        return setFloatValues(vector);
    }
    /** set/get an float values in this binary file */
    public boolean setLongValues(long[] vector)
    {
        if (dos==null)
        {
            StsException.systemError("StsBinaryFile.setLongValues() failed."
                + "File " + filename + " not properly opened for writing");
            return false;
        }

        try
        {
            if(includeSizeInFile)
                dos.writeInt(vector.length);  // save the size
            for (int i=0; i<vector.length; i++) dos.writeLong(vector[i]);
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.setLongValues() failed." +
                ": unable to write vector to file " + filename, e, StsException.WARNING);
            return false;
        }
    }

    /** set/get a float vector in this binary file */
    public boolean setLongVector(StsLongVector longVector)
    {
        if (dos==null)
        {
            StsException.systemError("StsBinaryFile.setLongVector() failed."
                + "File " + filename + " not properly opened for writing");
            return false;
        }

        try
        {
            if(longVector == null || longVector.getValues() == null)
            {
                dos.writeInt(0);
                dos.writeLong(0);
                dos.writeLong(0);
                return true;
            }

            long[] values = longVector.getValues();
            int nValues = values.length;
            dos.writeInt(nValues);  // save the size
            dos.writeLong(longVector.getMinValue());
            dos.writeLong(longVector.getMaxValue());

            for (int n=0; n<nValues; n++)
            {
                dos.writeLong(values[n]);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.setLongVector() failed." +
                "Unable to write vector to file " + filename, e, StsException.WARNING);
            return false;
        }

    }
    /** set/get an float values in this binary file */
    public boolean setLongValues(long[] vector, byte[] nullFlags, byte notNullFlag, long nullValue)
    {
        if (dos==null)
        {
            StsException.systemError("StsBinaryFile.setLongValues() failed."
                + "File " + filename + " not properly opened for writing");
            return false;
        }

        try
        {
            dos.writeInt(vector.length);  // save the size
            for (int i=0; i<vector.length; i++)
            {
                if(nullFlags[i] == notNullFlag)
                    dos.writeLong(vector[i]);
                else
                    dos.writeLong(nullValue);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.setLongValues() failed." +
                ": unable to write vector to file " + filename, e, StsException.WARNING);
            return false;
        }
    }
    /** set/get an float values in this binary file */
    public boolean setFloatValues(float[] vector)
    {
        if (dos==null)
        {
            StsException.systemError("StsBinaryFile.setFloatValues() failed."
                + "File " + filename + " not properly opened for writing");
			return false;
        }

        try
        {
            if(includeSizeInFile)
                dos.writeInt(vector.length);  // save the size
            for (int i=0; i<vector.length; i++) dos.writeFloat(vector[i]);
			return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.setFloatValues() failed." +
                ": unable to write vector to file " + filename, e, StsException.WARNING);
			return false;
        }
    }

    /** set/get a float vector in this binary file */
    public boolean setFloatVector(StsFloatVector floatVector)
    {
        if (dos==null)
        {
            StsException.systemError("StsBinaryFile.setFloatVector() failed."
                + "File " + filename + " not properly opened for writing");
			return false;
        }

        try
        {
            if(floatVector == null || floatVector.getValues() == null)
            {
                dos.writeInt(0);
                dos.writeFloat(0.0f);
                dos.writeFloat(0.0f);
                return true;
            }

			float[] values = floatVector.getValues();
			int nValues = values.length;
            dos.writeInt(nValues);  // save the size
			dos.writeFloat(floatVector.getMinValue());
			dos.writeFloat(floatVector.getMaxValue());
            for (int n=0; n<nValues; n++) dos.writeFloat(values[n]);
			return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.setFloatVector() failed." +
                "Unable to write vector to file " + filename, e, StsException.WARNING);
			return false;
        }

    }
    /** set/get an float values in this binary file */
    public boolean setFloatValues(float[] vector, byte[] nullFlags, byte notNullFlag, float nullValue)
    {
        if (dos==null)
        {
            StsException.systemError("StsBinaryFile.setFloatValues() failed."
                + "File " + filename + " not properly opened for writing");
			return false;
        }

        try
        {
            dos.writeInt(vector.length);  // save the size
            for (int i=0; i<vector.length; i++)
            {
                if(nullFlags[i] == notNullFlag)
                    dos.writeFloat(vector[i]);
                else
                    dos.writeFloat(nullValue);
            }
			return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.setFloatValues() failed." +
                ": unable to write vector to file " + filename, e, StsException.WARNING);
			return false;
        }
    }
    public float[] getFloatValues()
    {
        try
        {
            int nValues = dis.readInt(); // get the size
//            System.out.println("StsBinaryFile.getFloatValues() nValues: " + nValues);
            float[] values = new float[nValues];
            for (int i = 0; i < nValues; i++) values[i] = dis.readFloat();
            return values;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getFloatValues", e);
            return null;
        }
    }

    public boolean getFloatVector(StsFloatVector floatVector, boolean loadValues)
    {
        if (dis==null)
        {
            StsException.systemError("StsBinaryFile.getFloatVector() failed."
                + "File " + filename + " not properly opened for reading");
			return false;
        }

        try
        {
            int size = dis.readInt(); // get the size
			floatVector.setMinValue(dis.readFloat());
			floatVector.setMaxValue(dis.readFloat());

			if(loadValues)
			{
                float[] values = new float[size];
                for (int i=0; i<size; i++) values[i] = dis.readFloat();
                floatVector.setValues(values);
			}
			return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.getFloatVector() failed." +
                "Unable to read vector from file " + filename, e, StsException.WARNING);
			return false;
        }
    }
    public long[] getLongValues()
    {
        try
        {
            int nValues = dis.readInt(); // get the size
            long[] values = new long[nValues];
            for (int i = 0; i < nValues; i++) values[i] = dis.readLong();
            return values;
        }
        catch (Exception e)
        {
            StsException.outputException("getLongValues() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public boolean getLongVector(StsLongVector longVector, boolean loadValues)
    {
        if (dis==null)
        {
            StsException.systemError("StsBinaryFile.getLongVector() failed."
                + "File " + filename + " not properly opened for reading");
            return false;
        }

        try
        {
            int size = dis.readInt(); // get the size
            longVector.setMinValue(dis.readLong());
            longVector.setMaxValue(dis.readLong());
            if(loadValues)
            {
                long[] values = new long[size];
                for (int i=0; i<size; i++)
                {
                    values[i] = dis.readLong();
                }
                longVector.setValues(values);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.getLongVector() failed." +
                "Unable to read vector from file " + filename, e, StsException.WARNING);
            return false;
        }
    }
    /** set/get an float values in this binary file */
    public boolean setByteValues(byte[] vector, boolean incSize)
    {
        includeSizeInFile = incSize;
        return setByteValues(vector);
    }
    /** set/get a byte vector in this binary file */
    public boolean setByteValues(byte[] vector)
    {
        if (dos==null)
        {
            StsException.systemError("StsBinaryFile.setByteValues() failed."
                + "File " + filename + " not properly opened for writing");
			return false;
        }

        try
        {
            if(includeSizeInFile)
                dos.writeInt(vector.length);  // save the size
            for (int i=0; i<vector.length; i++) dos.writeByte(vector[i]);
			return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.setByteValues() failed." +
                ": unable to write vector to file " + filename, e, StsException.WARNING);
			return false;
        }
    }

    public byte[] getByteValues(int size)
    {
        int i = 0;

        if (dis==null)
        {
            StsException.systemError("StsBinaryFile.getByteValues() failed."
                + ": file " + filename + " not properly opened for reading");
			return null;
        }

        byte[] vector = null;
        try
        {
            vector = new byte[size];
            for (i=0; i<size; i++) vector[i] = dis.readByte();
            return vector;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.getByteValues() failed." +
                ": unable to read vector from file " + filename + ".\n" +
                i + " values read. Expected " + size,
                e, StsException.WARNING);
			return null;
        }
    }

    public byte[] getByteValues()
    {
        if (dis==null)
        {
            StsException.systemError("StsBinaryFile.getByteValues() failed."
                + ": file " + filename + " not properly opened for reading");
			return null;
        }

        byte[] vector = null;
        try
        {
            int size = dis.readInt();
            vector = new byte[size];
            for (int i=0; i<size; i++) vector[i] = dis.readByte();
            return vector;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.getByteValues() failed." +
                ": unable to read vector from file " + filename, e, StsException.WARNING);
			return null;
        }
    }

    /** set/get an integer vector in this binary file */
    public boolean setIntegerValues(int[] vector)
    {
        if (dos==null)
        {
            StsException.systemError("StsBinaryFile.setIntegerValues() failed."
                + "File " + filename + " not properly opened for writing");
			return false;
        }

        try
        {
            dos.writeInt(vector.length);  // save the size
            for (int i=0; i<vector.length; i++) dos.writeInt(vector[i]);
			return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.setIntegerValues() failed." +
                ": unable to write vector to file " + filename, e, StsException.WARNING);
			return false;
        }
    }

    public int[] getIntegerValues()
    {
        if (dis==null)
        {
            StsException.systemError("StsBinaryFile.getIntegerValues() failed."
                + ": file " + filename + " not properly opened for reading");
			return null;
        }

        int[] vector = null;
        try
        {
            int size = dis.readInt();
            vector = new int[size];
            for (int i=0; i<size; i++) vector[i] = dis.readInt();
            return vector;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.getIntegerValues() failed." +
                ": unable to read vector from file " + filename, e, StsException.WARNING);
			return null;
        }
    }

    /** set/get double values in this binary file */
    public boolean setDoubleValues(double[] vector)
    {
        if (dos==null)
        {
            StsException.systemError("StsBinaryFile.setDoubleValues() failed."
                + "File " + filename + " not properly opened for writing");
			return false;
        }

        try
        {
            dos.writeInt(vector.length);  // save the size
            for (int i=0; i<vector.length; i++) dos.writeDouble(vector[i]);
			return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.setDoubleValues() failed." +
                ": unable to write vector to file " + filename, e, StsException.WARNING);
			return false;
        }
    }

    public double[] getDoubleValues()
    {
        if (dis==null)
        {
            StsException.systemError("StsBinaryFile.getDoubleValues() failed."
                + ": file " + filename + " not properly opened for reading");
			return null;
        }

        double[] vector = null;
        try
        {
            int size = dis.readInt();
            vector = new double[size];
            for (int i=0; i<size; i++) vector[i] = dis.readDouble();
            return vector;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.getDoubleValues() failed." +
                ": unable to read vector from file " + filename, e, StsException.WARNING);
			return null;
        }
    }

    /** set/get boolean values in this binary file */
    public boolean setBooleanValues(boolean[] vector)
    {
        if (dos==null)
        {
            StsException.systemError("StsBinaryFile.setBooleanValues() failed."
                + "File " + filename + " not properly opened for writing");
			return false;
        }

        try
        {
            dos.writeInt(vector.length);  // save the size
            for (int i=0; i<vector.length; i++) dos.writeBoolean(vector[i]);
			return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.setBooleanValues() failed." +
                ": unable to write vector to file " + filename, e, StsException.WARNING);
			return false;
        }
    }

    public boolean[] getBooleanValues()
    {
        if (dis==null)
        {
            StsException.systemError("StsBinaryFile.getBooleanValues() failed."
                + ": file " + filename + " not properly opened for reading");
			return null;
        }

        boolean[] vector = null;
        try
        {
            int size = dis.readInt();
            vector = new boolean[size];
            for (int i=0; i<size; i++) vector[i] = dis.readBoolean();
            return vector;
        }
        catch (Exception e)
        {
            StsException.outputException("StsBinaryFile.getBooleanValues() failed." +
                ": unable to read vector from file " + filename, e, StsException.WARNING);
			return null;
        }
    }

   /** set/get boolean values in this binary file */
    public boolean setByteValue(byte value)
    {
        if (dos==null)
        {
            StsException.systemError(this, "setByteValue",
                "file " + filename + " not properly opened for writing");
			return false;
        }

        try
        {
            dos.writeByte(value);  // save the size
			return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "setByteValue",
                "unable to read vector from file " + filename, e);
			return false;
        }
    }

    public byte getByteValue()
    {
        if (dis==null)
        {
            StsException.systemError(this, "getByteValue",
                "file " + filename + " not properly opened for reading");
			return -1;
        }
        try
        {
            return dis.readByte();
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getByteValue",
                "unable to read vector from file " + filename, e);
			return -1;
        }
    }
}
