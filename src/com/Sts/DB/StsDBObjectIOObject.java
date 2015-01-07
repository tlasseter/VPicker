package com.Sts.DB;

import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.io.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2007</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsDBObjectIOObject implements StsDBObjectIO
{
	static final public boolean debug = Main.isDbCmdDebug;

	private StsDBTypeClass classType;

	public StsDBObjectIOObject(StsDBTypeObject classType)
	{
		this.classType = (StsDBTypeClass)classType;
	}

	/**
	  * copyObject
	  *
	  * @param oldObject Object
	  * @return Object
	  */
	 public Object copyObject(Object oldObject) throws IOException, IllegalAccessException
	 {
		 Object newObject = classType.newInstance();
		 if (newObject == null)
			 return null;
		 boolean skipIndexField = false;
		 classType.copyFields(oldObject, newObject, skipIndexField);
		 return newObject;
	 }

	 /**
	  * readObject
	  *
	  * @param in StsDBInputStream
	  * @return Object
	  * @throws IOException
	  */
	 public Object readObject(StsDBInputStream in) throws IOException
	 {
		 return readObject(in, null);
	 }

	 /**
	  * readObject
	  *
	  * @param in StsDBInputStream
	  * @param object Object
	  * @return Object
	  * @throws IOException
	  */
	 public Object readObject(StsDBInputStream in, Object object) throws IOException
	 {
		 try
		 {
			 if (object == null)
			 {
				 object = classType.newInstance();
				 if (debug)
					 System.out.println("StsDBObjectIOObject.read(DBInputStream, Object) object created. Class: Object " + StsToolkit.getSimpleClassname(object));
			 }
			 if (debug)
			 {
				 System.out.println("StsDBObjectIOObject::readObject reading object " + StsToolkit.getSimpleClassname(object));
             }

			 if (object == null)
				 return null;
			 StsDBField[] dbFields = classType.getDBFields();
			 if (dbFields == null)
				 return object;
			 int nFields = dbFields.length;
			 for (int i = 0; i < nFields; i++)
			 {
				 StsDBField dbField = dbFields[i];
				 if (debug)
				 {
					 System.out.println("   StsDBObjectIOObject::readObject reading field " + dbField.getFieldClassName() + "::" + dbField.getFieldName());
				 }
				 dbField.read(in, object);
			 }
			 return object;
		 }
		 catch (Exception e)
		 {
			 StsException.outputException("StsDBObjectIOObject.read(in, obj) failed.", e, StsException.WARNING);
			 return null;
		 }
	 }

	/**
	 * writeObject
	 *
	 * @param out StsDBOutputStream
	 * @param obj Object
	 * @throws IOException
	 * @throws IllegalAccessException
	 */
	public void writeObject(StsDBOutputStream out, Object obj) throws IOException, IllegalAccessException
	{
		if (debug)
		{
			StsException.systemDebug(this, "writeObject", "writing object " + StsToolkit.getSimpleClassname(obj));
		}
		StsDBField[] dbFields = classType.getDBFields();
		int nFields = dbFields.length;
		for (int i = 0; i < nFields; i++)
		{
			StsDBField dbField = dbFields[i];
			if (debug)
			{
				StsException.systemDebug(this, "writeObject", "writing field " + dbField.getFieldName() + " (" + dbField.getFieldClassName() + ")");
			}
			dbField.write(out, obj);
		}
	}

	public void exportObject(StsDBFileObjectTrader objectTrader, Object obj) throws IllegalAccessException
	{
		StsDBField[] dbFields = classType.getDBFields();
		int nFields = dbFields.length;
		for (int i = 0; i < nFields; i++)
		{
			StsDBField dbField = dbFields[i];
			dbField.export(objectTrader, obj);
		}
	}
}
