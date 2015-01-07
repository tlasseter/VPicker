
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DB.*;
import com.Sts.IO.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

/** This class defines String keys paired with integer values *
  * These need to be saved in the database as they define the state of
  * any flags of interest
  */

public class StsProperties implements Serializable, StsCustomSerializable
{
    transient Hashtable properties;
    public static final String Separator = "=";

	public StsProperties()
	{
        this(10);
	}

	public StsProperties(int initSize)
	{
        properties = new Hashtable(initSize);
	}

    /** Get the value associated with this key: set value to 0 if not found */
    public int getInteger(String name)
    {
        String value = null;

        try
        {
            value = (String)properties.get(name);
            if(value == null) return 0;
            return Integer.parseInt(value);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getInteger", "Returned: " + value, e);
            return 0;
        }
    }

    /** Get the value associated with this key: set value to 0 if not found */
    public byte getByte(String name)
    {
        String value = null;

        try
        {
            value = (String)properties.get(name);
            if(value == null) return 0;
            return Byte.parseByte(value);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getByte", "Returned: " + value, e);
            return 0;
        }
    }

    /** check if this proerty exists */
    public boolean hasProperty(String name)
    {
        return properties.get(name) != null;
    }

    /** return the value as a boolean */
    public boolean getBoolean(String name)
    {
        String value = null;

        try
        {
            value = (String)properties.get(name);
            if(value == null) return false;
            return value.equalsIgnoreCase("true");
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getBoolean", "Returned: " + value, e);
            return false;
        }
    }

    public String getProperty(String name)
    {
        return (String)properties.get(name);
    }

    /** Set the value associated with this key: if key doesn't exist add key and value */
    public void set(String name, int value)
    {
        properties.put(name, Integer.toString(value));
    }

    /** Set the key value given a boolean */
    public void set(String name, boolean value)
    {
        try
        {
    	    properties.put(name, new Boolean(value).toString());
        }
        catch(Exception e)
        {
            StsException.outputException("StsProperties.set(boolean) failed.",
                e, StsException.WARNING);
        }
    }

    public void set(String name, String value)
    {
		if(name == null || value == null) return;
        properties.put(name, value);
    }

    public String getKeyFromValue(String value)
    {
        Object[] values = properties.values().toArray();
        Object[] keys = properties.keySet().toArray();
        for(int n = 0; n < values.length; n++)
            if(values[n].equals(value)) return (String)keys[n];
        return null;
    }

    public String getKeyFromValue(int i)
    {
        return getKeyFromValue(Integer.toString(i));
    }

    public Object[] getKeysAsArray()
    {
        return (Object[])properties.keySet().toArray();
    }

    public Object[] getValuesAsArray()
    {
        return (Object[])properties.values().toArray();
    }

    public void print()
    {
        Enumeration enumeration = properties.keys();
        while(enumeration.hasMoreElements())
        {
            Object key = enumeration.nextElement();
            String value = properties.get(key).toString();
        	System.out.println("key: " + key + " value: " + value);
        }
    }

	public void writeObject(StsDBOutputStream out) throws IllegalAccessException, IOException
	{
  	    Object[] entries = properties.entrySet().toArray();
		int nEntries = entries.length;
        out.writeInt(nEntries);
        for(int n = 0; n < nEntries; n++)
		{
			Map.Entry entry = (Map.Entry)entries[n];
            out.writeUTF(entry.getKey().toString());
            out.writeUTF(entry.getValue().toString());
        }
	}

	public void readObject(StsDBInputStream in) throws IllegalAccessException, IOException
	{
        int nEntries = in.readInt();
        properties = new Hashtable(nEntries);
        for(int n = 0; n < nEntries; n++)
        {
            String key = in.readUTF();
            String value = in.readUTF();
            properties.put(key, value);
        }
	}

    public void exportObject(StsDBFileObjectTrader objectTrader) { } 
       public static void main(String[] args)
    {
        StsProperties properties = new StsProperties(10);

        String booleanKey = new String("Boolean true:");
        properties.set(booleanKey, true);
        String stringKey = new String("aaaa");
        properties.set(stringKey, "aaaa");
        String intKey = new String("Int 5:");
        properties.set(intKey, 5);

        properties.print();

        System.out.println(booleanKey + " " + properties.getBoolean(booleanKey));
        System.out.println(stringKey + " " + properties.getProperty(stringKey));
        System.out.println(intKey + " " + properties.getInteger(intKey));

	    Object[] entries = properties.properties.entrySet().toArray();
		int nEntries = entries.length;
		for(int n = 0; n < nEntries; n++)
		{
			Map.Entry entry = (Map.Entry)entries[n];
			System.out.println("KEY: " + entry.getKey().toString() + " VALUE: " + entry.getValue().toString());
		}
    }

	public void writeToAsciiFile(StsAsciiFile parmFile) throws IOException {
		parmFile.openWrite();
		Object[] entries = properties.entrySet().toArray();
		for (int i=0; i< entries.length; i++) {
			Map.Entry entry = (Map.Entry)entries[i];
			parmFile.writeLine(entry.getKey()+Separator+entry.getValue());
		}
		parmFile.close();
	}

	public void readFromFile(StsAsciiFile parmFile) throws IOException {
		parmFile.openRead();
		String line = parmFile.readLine();
		while (line != null) {
			String[] split = line.split(Separator);
			if (split.length == 2) {
				String key = split[0].trim();
				String val = split[1].trim();
				set(key, val);
			}
			line = parmFile.readLine();
		}
		parmFile.close();
	}
}

