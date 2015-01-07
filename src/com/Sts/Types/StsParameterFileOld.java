package com.Sts.Types;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.lang.reflect.*;
import java.util.*;

public class StsParameterFileOld
{
    static protected ArrayList missingFields = null;
	static private final String stringClass = "java.lang.String";

	static private final boolean debug = false;

    static public boolean readObjectFields(String pathname, Object obj) throws StsException
    {
        Class c = obj.getClass();
        return readObjectFields(pathname, obj, c, c);
	}

	static public void writeObjectFields(String pathname, Object obj) throws StsException
	{
		Class c = obj.getClass();
		writeObjectFields(pathname, obj, c, c);
	}

  	/** write fields packed in a string to a file */
  	static public void writeObjectFields(String pathname, Object obj, Class subClass, Class superClass) throws StsException
  	{
        StsAsciiFile asciiFile = null;

    	if(pathname == null)
        {
      		throw new StsException(StsException.WARNING, "StsFieldIO.writeObjectFields:  Pathname is null");
    	}

    	try
        {
            StsFile file = StsFile.constructor(pathname);
			if(file.exists()) file.delete();
            asciiFile = new StsAsciiFile(file);
            asciiFile.openWrite();
            writeObjectFields(asciiFile, subClass, superClass, obj, "");
    	}
    	catch (Exception ex)
        {
        	StsException.outputException("StsFieldIO.writeObjectFields: ", ex, StsException.WARNING);

    	}
    	finally
        {
      		try { if(asciiFile != null) asciiFile.close(); }
      		catch(Exception ex) { }
    	}
  	}
    /** output field name-value pairs from defined class up thru and including the superclass.
     *  format the object fields as a strings
     */
    static public void writeObjectFields(StsAsciiFile file, Class subClass, Class superClass, Object object, String prefix) throws StsException
    {
        try
        {
            if(subClass == null)
                subClass = object.getClass();

            Object[] fields = StsToolkit.getFields(object, subClass, superClass, false);
            int nFields = fields.length;
            for(int n = 0; n < nFields; n++)
            {
                Field field = (Field)fields[n];
                if(field.get(object) == null) continue;
				if(Modifier.isTransient(field.getModifiers())) continue;

                String fieldName = field.getName();
                Class type = field.getType();
                String typeName = type.getName();
                Object fieldObject = field.get(object);
                if(type.isPrimitive() || type == String.class)
                    file.writeLine(prefix + fieldName + " " + fieldObject.toString());
                else if(type.isArray() )
                {
                    file.writeString(prefix + fieldName);

                    Class arrayType = type.getComponentType();
					if(arrayType.isArray())
					{
						StsException.systemError("Developer!  StsParameterFileOld cannot handle multidimensioned arrays for field " + fieldName + " in object " + StsToolkit.getSimpleClassname(object) + "!  Sorry!");
						return;
					}
                    if(arrayType.equals(Integer.TYPE))
                    {
                        int[] ints = (int[])fieldObject;
                        for(int i = 0; i < ints.length; i++)
                            file.writeString(Integer.toString(ints[i]));
                    }
                    else if(arrayType.equals(Float.TYPE))
                    {
                        float[] floats = (float[])fieldObject;
                        for(int i = 0; i < floats.length; i++)
                            file.writeString(Float.toString(floats[i]));
                    }
					else if(arrayType.equals(Double.TYPE))
					{
						double[] doubles = (double[])fieldObject;
						for(int i = 0; i < doubles.length; i++)
							file.writeString(Double.toString(doubles[i]));
				    }
                    else if(arrayType.equals(Boolean.TYPE))
                    {
                        boolean[] booleans = (boolean[])fieldObject;
                        for(int i = 0; i < booleans.length; i++)
                            file.writeString(Boolean.toString(booleans[i]));
				    }
					else if(arrayType.equals(Byte.TYPE))
					{
						byte[] bytes = (byte[])fieldObject;
						for(int i = 0; i < bytes.length; i++)
							file.writeString(Byte.toString(bytes[i]));
				    }
					else if(arrayType.equals(Character.TYPE))
					{
						char[] chars = (char[])fieldObject;
						for(int i = 0; i < chars.length; i++)
							file.writeString(Character.toString(chars[i]));
				    }
					else if(arrayType.equals(Short.TYPE))
					{
						short[] shorts = (short[])fieldObject;
						for(int i = 0; i < shorts.length; i++)
							file.writeString(Short.toString(shorts[i]));
				    }
					else if(arrayType.equals(Long.TYPE))
					{
						long[] longs = (long[])fieldObject;
						for(int i = 0; i < longs.length; i++)
							file.writeString(Long.toString(longs[i]));
				    }
					else if(arrayType.equals(String.class))
					{
						  String[] strings = (String[])fieldObject;
						   for(int i = 0; i < strings.length; i++)
						   {
							   strings[i].replaceAll(" ", "");
							   file.writeString(strings[i]);
						   }
                    }
					else
					{
						if (fieldObject != null && fieldObject instanceof Object[])
						{
							Object[] objects = (Object[]) fieldObject;
							int nObjects = objects.length;
							file.writeLine(Integer.toString(nObjects));
							for (int i = 0; i < nObjects; i++)
								writeObjectFields(file, type, type, objects[i], prefix + fieldName + ".");
						}
					}
                    file.endLine();
                }
                else
                {
                    writeObjectFields(file, null, null, field.get(object), prefix + fieldName + ".");
                }
            }
        }
        catch(Exception e)
        {
            throw new StsException(StsException.WARNING, e.toString());
        }
  	}

    /** fill the object fields with values from a string */
    static public boolean unpackObjectFields(Object object, Class superClass, String s) throws StsException
    {
        boolean isMissingFields = false;

        Class subClass = object.getClass();
        Object[] fields = StsToolkit.getFields(object, subClass, superClass, false);
        int nFields = fields.length;
        HashMap fieldsHashMap = new HashMap(nFields);

        for(int n = 0; n < nFields; n++)
        {
            Field field = (Field)fields[n];
            fieldsHashMap.put(field.getName(), field);
        }

        missingFields = new ArrayList();
        StringTokenizer stok = new StringTokenizer(s);
        int nTokenPairs = stok.countTokens() / 2;

        for (int i=0; i<nTokenPairs; i++)
        {
			String fieldName = stok.nextToken();
        	String fieldVal = stok.nextToken();
            Field field = (Field)fieldsHashMap.get(fieldName);
            if(field != null)
            {
                setObjectFieldValue(object, field, fieldVal, fieldName);
        	}
            else
            {
                isMissingFields = true;
                missingFields.add(fieldName);
//          		throw new StsException(StsException.WARNING, "StsProjParms.unpackFields:" +
//                        " Missing value for required parameter '" + fields[i].getName() + "'");
            }
      	}
        return isMissingFields;
  	}

    /** fill the object fields with values from a string */
    static public boolean readObjectFields(String pathname, Object object, Class subClass, Class superClass) throws StsException
    {
        StsAsciiFile asciiFile = null;

        try
        {
            StsFile file = StsFile.constructor(pathname);
            if(!file.exists()) return false;
            asciiFile = new StsAsciiFile(file);
            if(!asciiFile.openReadWithErrorMessage()) return false;
            return readObjectFields(asciiFile, object, subClass, superClass);
        }
        catch(Exception e)
        {
            StsException.outputException("StsParameterFileOld.readObjectFields() failed.",
                e, StsException.WARNING);
            return false;
        }
        finally
        {
      		try { if(asciiFile != null) asciiFile.close(); }
      		catch(Exception ex) { }
        }
    }

    /** fill the object fields with values from a string
     *  return false if any fields are missing
     */
    static public boolean readObjectFields(StsAsciiFile file, Object object, Class subClass, Class superClass) throws IOException, StsException
    {
        boolean isMissingFields = false;
        TreeMap fieldsTreeMap = StsToolkit.getFieldsTreeMap(object, subClass, superClass, false);
        missingFields = new ArrayList();
        String line;
        String fieldName, fieldVal;

        while((line = file.readLine()) != null)
        {
            StringTokenizer stok = new StringTokenizer(line);
            int nTokenPairs = stok.countTokens();
            if(nTokenPairs < 2) continue;

            fieldName = stok.nextToken();
            if(nTokenPairs == 2)
            {
                fieldVal = stok.nextToken();
            }
            else
            {
                String secondToken = stok.nextToken();
                int startIndex = line.indexOf(secondToken);
                int endIndex = line.length();
                fieldVal = line.substring(startIndex, endIndex);
            }
            Field field = (Field)fieldsTreeMap.get(fieldName);
            if(field != null)
            {
                setObjectFieldValue(object, field, fieldVal, fieldName);
            }
            else
            {
                isMissingFields = true;
                missingFields.add(fieldName);
//                System.out.println("Parameter file " + file.createFilename() + " missing field " + fieldName);
//          		throw new StsException(StsException.WARNING, "StsProjParms.unpackFields:" +
//                        " Missing value for required parameter '" + fields[i].getName() + "'");
            }
        }
        return !isMissingFields;
  	}

    /** format the object fields as a string */
    public String packFields() throws StsException
    {
        try
        {
            StringBuffer sb = new StringBuffer();
            Field[] fields = getClass().getDeclaredFields();
            for (int i = 0; i < fields.length; i++)
            {
                if (fields[i].get(this)==null ||
                    Modifier.isTransient(fields[i].getModifiers())) continue;
                sb.append(fields[i].getName() + " " +
                    fields[i].get(this).toString() + "\n");
            }
            return sb.toString();
        }
        catch(Exception e)
        {
            throw new StsException(StsException.WARNING, e.toString());
        }
  	}

    /** fill the object fields with values from a string */
    public void unpackFields(String s) throws StsException
    {
        Class thisClass = getClass();
        Field[] fields = thisClass.getDeclaredFields();
        int nFields = fields.length;
        boolean[] fieldFound = new boolean[nFields];
        for (int i=0; i<nFields; i++) fieldFound[i] = false;

        StringTokenizer stok = new StringTokenizer(s);
        int nTokenPairs = stok.countTokens() / 2;
        String fieldName, fieldVal;
        Field field;

        for (int i=0; i<nTokenPairs; i++)
        {
			fieldName = stok.nextToken();
        	fieldVal = stok.nextToken();
        	for (int j=0; j<nFields; j++)
            {
          		if (fieldName.equals(fields[j].getName()))
                {
            		setFieldValue(fields[j], fieldVal);
            		fieldFound[j] = true;
                    break;
          		}
        	}
      	}
      	for (int i=0; i<nFields; i++)
        {
            // see if the field was found and if it's transient (optional)
        	if (!fieldFound[i] && !Modifier.isTransient(fields[i].getModifiers()))
            {
          		throw new StsException(StsException.FATAL, "StsProjParms.unpackFields:" +
                        " Missing value for required parameter '" + fields[i].getName() + "'");

        	}
      	}
  	}

  	static private void setObjectFieldValue(Object object, Field f, String val, String fieldName) throws StsException
  	{
    	try
        {
            f.setAccessible(true);

      		Class type = f.getType();
			if(debug) System.out.println("setting field " + fieldName + " with value " + val + " of type " + type.getName());
            if(type.isPrimitive())
            {
                if (type.equals(Double.TYPE)) f.setDouble(object, Double.valueOf(val).doubleValue());
                else if (type.equals(Float.TYPE)) f.set(object, (Object)Float.valueOf(val));
                else if (type.equals(Integer.TYPE)) f.setInt(object, Integer.valueOf(val).intValue());
                else if (type.equals(Long.TYPE)) f.setLong(object, Long.valueOf(val).longValue());
                else if (type.equals(Short.TYPE)) f.setShort(object, Short.valueOf(val).shortValue());
                else if (type.equals(Boolean.TYPE)) f.setBoolean(object, Boolean.valueOf(val).booleanValue());
                else if (type.equals(Byte.TYPE)) f.setByte(object, Byte.valueOf(val).byteValue());
                else if (type.equals(Character.TYPE)) f.setChar(object, val.charAt(0));
            }
			else if (type.getName().equals("java.lang.String")) f.set(object, (Object)val);
            else if(type.isArray())
            {
                StringTokenizer tokens = new StringTokenizer(val);
                int nTokens = tokens.countTokens();
                type = type.getComponentType();
                String arrayName = type.getName();
                int n = 0;
				 if(type.isArray())
				 {
					 StsException.systemError("Developer!  StsParameterFileOld cannot handle multidimensioned arrays for field " + fieldName + " in object " + StsToolkit.getSimpleClassname(object) + "!  Sorry!");
					 return;
				 }
				 if(type.equals(Integer.TYPE))
				 {
					 int[] ints = new int[nTokens];
					 while(tokens.hasMoreElements())
						 ints[n++] = Integer.parseInt(tokens.nextToken());
					 f.set(object, ints);
				 }
				 else if(type.equals(Float.TYPE))
				 {
					 float[] floats = new float[nTokens];
					 while(tokens.hasMoreElements())
						 floats[n++] = Float.parseFloat(tokens.nextToken());
					 f.set(object, floats);
				 }
				 else if(type.equals(Double.TYPE))
				 {
					 double[] doubles = new double[nTokens];
					  while(tokens.hasMoreElements())
						  doubles[n++] = Double.parseDouble(tokens.nextToken());
					f.set(object, doubles);
				 }
				 else if(type.equals(Boolean.TYPE))
				 {
					 boolean[] booleans = new boolean[nTokens];
					 while(tokens.hasMoreElements())
						 booleans[n++] = Boolean.valueOf(tokens.nextToken()).booleanValue();
					 f.set(object, booleans);
				 }
				 else if(type.equals(Byte.TYPE))
				 {
					 byte[] bytes = new byte[nTokens];
					 while(tokens.hasMoreElements())
						 bytes[n++] = Byte.parseByte(tokens.nextToken());
					 f.set(object, bytes);
				 }
				 else if(type.equals(Character.TYPE))
				 {
					 char[] chars = new char[nTokens];
					 while(tokens.hasMoreElements())
						 chars[n++] = tokens.nextToken().charAt(0);
					 f.set(object, chars);
				 }
				 else if(type.equals(Short.TYPE))
				 {
					 short[] shorts = new short[nTokens];
					 while(tokens.hasMoreElements())
						 shorts[n++] = Short.parseShort(tokens.nextToken());
					 f.set(object, shorts);
				 }
				 else if(type.equals(Long.TYPE))
				 {
					 long[] longs = new long[nTokens];
					 while(tokens.hasMoreElements())
						 longs[n++] = Long.parseLong(tokens.nextToken());
					 f.set(object, longs);
				 }
				 else if(type.equals(String.class))
				 {
					 String[] strings = new String[nTokens];
					 while(tokens.hasMoreElements())
						 strings[n++] = tokens.nextToken();
					 f.set(object, strings);
				 }
				 else // can't determine type; treat as Strings
				 {
					 StsException.systemError("StsParameterFileOld.setObjectFieldValue() failed. Couldn't fine class type: " + type.getName());
				 }
            }
    	}
    	catch (IllegalAccessException e)
        {
      		throw new StsException(StsException.WARNING, "StsProjParms.setFieldValue: " + e.toString() +
                                    "\nFor field: " + fieldName);
    	}
  	}
/*
  	static private void setObjectFieldValue(Object object, Field f, String val, String fieldName) throws StsException
  	{
    	try
        {
            f.setAccessible(true);

      		Class type = f.getType();
			boolean isaString = type.equals(String.class);
			if(mainDebug) System.out.println("setting field " + fieldName + " with value " + val + " of type " + type.getName());
            if(type.isPrimitive())
            {
                if (type.equals(Double.TYPE)) f.setDouble(object, Double.valueOf(val).doubleValue());
                else if (type.equals(Float.TYPE)) f.set(object, (Object)Float.valueOf(val));
                else if (type.equals(Integer.TYPE)) f.setInt(object, Integer.valueOf(val).intValue());
                else if (type.equals(Long.TYPE)) f.setLong(object, Long.valueOf(val).longValue());
                else if (type.equals(Short.TYPE)) f.setShort(object, Short.valueOf(val).shortValue());
                else if (type.equals(Boolean.TYPE)) f.setBoolean(object, Boolean.valueOf(val).booleanValue());
                else if (type.equals(Byte.TYPE)) f.setByte(object, Byte.valueOf(val).byteValue());
                else if (type.equals(Character.TYPE)) f.setChar(object, val.charAt(0));
            }
			else if (type.getName().equals("java.lang.String"))
				f.set(object, (Object)val);
//            else if (type.equals(String.class)) f.set(object, val);
            else if(type.isArray())
            {
                StringTokenizer tokens = new StringTokenizer(val);
                int nTokens = tokens.countTokens();
                type = type.getComponentType();
				if(type.isArray())
				{
					StsException.systemError("Developer!  StsParameterFileOld cannot handle multidimensioned arrays!  Sorry!");
					return;
				}
                int n = 0;
                if(type.equals(Integer.TYPE))
                {
                    int[] ints = new int[nTokens];
                    while(tokens.hasMoreElements())
                        ints[n++] = Integer.parseInt(tokens.nextToken());
                    f.set(object, ints);
                }
                else if(type.equals(Float.TYPE))
                {
                    float[] floats = new float[nTokens];
                    while(tokens.hasMoreElements())
                        floats[n++] = Float.parseFloat(tokens.nextToken());
                    f.set(object, floats);
                }
				else if(type.equals(Double.TYPE))
				{
					double[] doubles = new double[nTokens];
					while(tokens.hasMoreElements())
						doubles[n++] = Double.parseDouble(tokens.nextToken());
					f.set(object, doubles);
				}
                else if(type.equals(Boolean.TYPE))
                {
                    boolean[] booleans = new boolean[nTokens];
                    while(tokens.hasMoreElements())
                        booleans[n++] = Boolean.valueOf(tokens.nextToken()).booleanValue();
                    f.set(object, booleans);
                }
				else if(type == String.class)
				{
					String[] strings = new String[nTokens];
					while(tokens.hasMoreElements())
						strings[n++] = tokens.nextToken();
                        f.set(object, strings);
                }
				else if(type.equals(Character.TYPE))
				{
					 char[] chars = new char[nTokens];
					 while(tokens.hasMoreElements())
						 chars[n++] = tokens.nextToken().charAt(0);
					 f.set(object, chars);
                }
				else if(type.equals(Short.TYPE))
				{
					short[] shorts = new short[nTokens];
					while(tokens.hasMoreElements())
						shorts[n++] = Short.parseShort(tokens.nextToken());
					f.set(object, shorts);
				}
				else if(type.equals(Long.TYPE))
				{
					long[] longs = new long[nTokens];
					while(tokens.hasMoreElements())
						longs[n++] = Long.parseLong(tokens.nextToken());
					f.set(object, longs);
				}
            }
    	}
    	catch (Exception e)
        {
      		throw new StsException(StsException.WARNING, "StsProjParms.setFieldValue: " + e.toString() +
                                    "\nFor field: " + fieldName);
    	}
  	}
*/
    /* parse a value from a string and set it in the right field */
  	private void setFieldValue(Field f, String val) throws StsException
  	{
    	try
        {
      		Class type = f.getType();
      		if (type.equals(Double.TYPE)) {
        		f.setDouble(this, Double.valueOf(val).doubleValue());
            }
            else if (type.equals(Float.TYPE)) {
            	f.set(this, (Object)Float.valueOf(val));
            }
            else if (type.equals(Integer.TYPE)) {
            	f.setInt(this, Integer.valueOf(val).intValue());
            }
            else if (type.equals(Long.TYPE)) {
            	f.setLong(this, Long.valueOf(val).longValue());
            }
            else if (type.equals(Short.TYPE)) {
            	f.setShort(this, Short.valueOf(val).shortValue());
            }
            else if (type.equals(Boolean.TYPE)) {
            	f.setBoolean(this, Boolean.valueOf(val).booleanValue());
            }
            else if (type.equals(Byte.TYPE)) {
            	f.setByte(this, Byte.valueOf(val).byteValue());
            }
            else if (type.equals(Character.TYPE)) {
            	f.setChar(this, val.charAt(0));
            }
            else if (type.getName().equals("java.lang.String")) {
            	f.set(this, (Object)val);
            }
    	}
    	catch (IllegalAccessException e)
        {
      		throw new StsException(StsException.FATAL, "StsProjParms.setFieldValue: " + e.toString());
    	}
  	}
  	/** read serialized parameters from a file into a new object */
    static public Object readSerializedFields(String filename) throws StsException
  	{
    	FileInputStream fileIS = null;
    	try
        {
      		fileIS = new FileInputStream(filename);
      		ObjectInputStream objIS= new ObjectInputStream(fileIS);
      		return objIS.readObject();
    	}
    	catch (Exception ex)
        {
//      	throw new StsException(StsException.FATAL, "StsFieldIO.readSerializedFields: " + ex.toString());
			StsException.outputException("StsFieldIO.readSerializedFields: ", ex, StsException.WARNING);
            return null;
     	}
    	finally
        {
      		try
            {
        		fileIS.close();
      		}
      		catch (IOException ex)
        	{
//      		throw new StsException(StsException.FATAL, "StsFieldIO.readSerializedFields: " + ex.toString());
         		StsException.outputException("StsFieldIO.readSerializedFields: ", ex, StsException.FATAL);
				return null;
      		}
    	}
  	}

  	/** write serialized parameters into a file */
  	static public void writeSerializedFields(Object obj, String filename) throws StsException
  	{
    	FileOutputStream fileOS = null;
    	try
        {
      		fileOS = new FileOutputStream(filename);
      		ObjectOutput objOut = new ObjectOutputStream(fileOS);
      		objOut.writeObject(obj);
      		objOut.flush();
    	}
    	catch (Exception ex)
        {
        	StsException.outputException("StsFieldIO.writeSerializedFields: ", ex, StsException.FATAL);
    	}
    	finally
        {
      		try
            {
        		fileOS.close();
      		}
      		catch (IOException ex)
            {
        		StsException.outputException("StsFieldIO.writeSerializedFields: ", ex, StsException.FATAL);
      		}
    	}
  	}

  	/** read serialized fields into to a new object with filename prompt */
  	static public Object readSerializedFields() throws StsException
  	{
    	String filename = selectFileLoad(null, null);
    	if (filename != null)  // If user didn't click "Cancel."
        {
      		return readSerializedFields(filename);
    	}
    	return (Object)null;
  	}

  	/** write serialized fields from an object with filename prompt */
  	static public boolean writeSerializedFields(Object obj) throws StsException
  	{
    	String filename = selectFileSave(null);
    	if (filename != null)  // If user didn't click "Cancel."
        {
      		writeSerializedFields(obj, filename);
      		return true;
    	}
    	return false;
  	}

  	/** get a file name */
  	static public File getSaveFile(String title, String directory, String filter, String defaultName)
  	{
    	// Create a file dialog to query the user for a filename.
        StsWin3d win3d = StsObject.getCurrentModel().win3d;
    	Frame frame = win3d != null ? (Frame)win3d : new Frame();
        StsFileChooser chooser = StsFileChooser.createFileChooserPrefix(frame, title, directory, filter);
//		StsFileChooser chooser = new StsFileChooser(frame, title, directory, filter, false, false);
		chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
		chooser.getFileChooser().setName(defaultName);
        while(chooser.show())
        {
            File f = new File(chooser.getFilePath());
            if(!f.getName().endsWith(filter))
            	f = new File(chooser.getFilePath() + ".csv");

            if(!f.exists())
                return f;
            else
            {
            	if(StsYesNoDialog.questionValue(frame, "File exist, do you want to overwrite it?"))
            		return f;
            }
        }
        return null;
  	}

  	/** get a file name */
  	static public String selectFile(String title, int mode, String directory, String filter, String filterString)
  	{
    	// Create a file dialog to query the user for a filename.
        StsWin3d win3d = StsObject.getCurrentModel().win3d;
    	Frame frame = win3d != null ? (Frame)win3d : new Frame();
        StsFileChooser chooser = StsFileChooser.createFileChooserPrefix(frame, title, directory, filter);
		//StsFileChooser chooser = new StsFileChooser(frame, title, directory, filter, true, false);
        String dir = null;
        while( chooser.show() )
        {
            File f = new File(chooser.getFilePath());
            if( !f.exists() )
                StsMessageFiles.logMessage("File " + f.getPath() + " not found.");
            else if( f.isDirectory() )
                StsMessageFiles.logMessage("Selected file " + f.getPath() + " is a directory.");
            else return f.getPath();
        }
        return null;
  	}

  	/** get a file name */
  	static public String selectFileSave(String dir)
  	{
    	return selectFile("Select a file to save", FileDialog.SAVE, dir, null, null);
  	}

  	/** get a file name */
  	static public String selectFileLoad(String dir, String filter)
  	{
    	return selectFile("Select a file to open", FileDialog.LOAD, dir, filter, " ");
  	}

  	static public String selectFileLoad(String dir, String filter, String filterString)
  	{
    	return selectFile("Select a file to open", FileDialog.LOAD, dir, filter, filterString);
  	}

  	/** read fields packed in a string into an object with filename prompt */
  	static public String askReadPackedFields(Object obj, String filter, String filterString) throws StsException
  	{
    	String filename = null;
        while( true )
        {
    	    filename = selectFileLoad(null, filter, filterString);
            if( filename == null ) break;
            File f = new File(filename);
//            if( !StsProjectParms.isProjectFile(f) )
//                StsMessageFiles.logMessage("Selected file " + f.getPath() + " is not a valid S2S project file.");
//            else break;
        }
    	if (filename != null)  // If user didn't click "Cancel."
        {
      		readObjectFields(filename, obj, null, null);
      		return filename;
    	}
    	return null;
  	}

  	static public String askReadPackedFields(Object obj, String filter) throws StsException
  	{
    	String filename = selectFileLoad(null, filter);
    	if (filename != null)  // If user didn't click "Cancel."
        {
      		readObjectFields(filename, obj, null, null);
      		return filename;
    	}
    	return null;
  	}

    static public void main(String[] args)
    {
        try
        {

			StsPreStackLine3d line = new StsPreStackLine3d();
			StsParameterFileOld.readObjectFields("j:/TsunamiSEGY/prestack3d.txt.summit_tape01.line.0", line, StsPreStackLine.class, StsBoundingBox.class);
			System.out.println("stemname: " + line.stemname);
		/*
		    TestSub writeTest = TestSub.writeTest();
			writeTest.print("write");
            StsParameterFileOld.writeObjectFields("c:/test.txt", writeTest, Test.class, TestSuper.class);
            TestSub readTest = TestSub.readTest();
            StsParameterFileOld.readObjectFields("c:/test.txt", readTest, Test.class, TestSuper.class);
            readTest.print("read");
        */
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }
}

class MemberOld
{
    int mi = 21;
    float mf = 22.0f;

    MemberOld()
    {
    }

    void initialize()
    {
        mi = 21;
        mf = 22.0f;
    }

    void print()
    {
        System.out.println("mi = " + mi + " mf = " + mf);
    }
}

abstract class TestSuperOld
{
	String superString;
	float superF;
}

class TestSubOld extends TestOld
{
	static TestSub writeTest()
	{
		TestSub test = new TestSub();
		test.initialize();
		return test;
	}

	static TestSubOld readTest()
	{
		return new TestSubOld();
    }
}

abstract class TestOld extends TestSuperOld
{
    String string;
    int i;
    float f;
    double d;
    byte b;
    float[] ff;
    MemberOld member;

    TestOld()
    {
		super();
    }

    void initialize()
    {
		string = "string";
		superString = "superString";
		superF = 99.0f;
        i = 1;
        f = 2.0f;
        d = 3;
        b = 4;
        ff = new float[] { 11f, 12f };
        member = new MemberOld();
        member.initialize();
    }

    void print(String printString)
    {
		System.out.println(printString);
        System.out.println("string: '" + string + "'");
		System.out.println("superString: '" + superString + "'");
		System.out.println("superF: " + superF);
        System.out.println("i: " + i);
        System.out.println("f: " + f);
        System.out.println("d: " + d);
        System.out.println("b: " + b);
        if(ff != null)
        {
            System.out.print("ff:");
            for(int n = 0; n < ff.length; n++)
                System.out.print(" " + ff[n]);
            System.out.println();
        }
        if(member != null) member.print();
    }
/*
    static Test writeTest()
    {
        Test test = new Test();
        test.classInitialize();
        return test;
    }

    static Test readTest()
    {
        return new Test();
    }
*/
}
