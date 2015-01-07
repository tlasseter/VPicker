package com.Sts.DB;

import com.Sts.DB.DBCommand.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.lang.reflect.*;

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

public class StsDBFileObject extends StsDBFile
{
	private StsDBFileObject(StsAbstractFile file, StsProgressBar progressBar) throws StsException, IOException
	{
		super(file, progressBar);
	}

	static public StsDBFileObject openWrite(StsAbstractFile file, StsProgressBar processView)
	{
		try
		{
			StsDBFileObject dbFile = new StsDBFileObject(file, processView);
			if (!dbFile.openWrite())
			   return null;
			return dbFile;
		}

		catch (Exception e)
		{
			new StsMessage(null, StsMessage.ERROR, "StsDBFile.openWrite() failed. File: " + file.getFilename());
			return null;
		}
	}

	static public StsDBFileObject openRead(StsAbstractFile file, StsProgressBar processView)
	{
		try
		{
			StsDBFileObject dbFile = new StsDBFileObject(file, processView);
			if (! dbFile.openReadAndCheckFile(true))
			{
				return null;
			}
			return dbFile;
		}

		catch (Exception e)
		{
			new StsMessage(null, StsMessage.ERROR, "StsDBFile.openReadAndCheck() failed. File: " + file.getFilename());
			return null;
		}
	}

	static public boolean writeObjectFile(String pathname, Object object, StsProgressBar processView)
	{
		try
		{
			StsFile file = StsFile.constructor(pathname);
			StsDBFileObject currentDBFile = StsDBFileObject.openWrite(file, processView);
//			currentDBFile.commitCmd("Reference Model Transaction", new StsModelDBRefCmd());
			currentDBFile.commitCmd("Save Object Transaction", new StsSimpleObjectAddCmd(object));
			currentDBFile.close();
			currentDBFile = null;
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsDBFile.writeObjectFile() failed.", e, StsException.WARNING);
			return false;
		}
	}

    // read an Object which is tied to this database; transactionTimes will be checked against DB to make sure they are in sync
    static public boolean readDatabaseObjectFile(String pathname, Object object, StsProgressBar processView)
	{
		try
		{
			StsFile file = StsFile.constructor(pathname);
			StsDBFileObject currentDBFile = StsDBFileObject.openRead(file, processView);
			StsSimpleObjectAddCmd.setObject(object);
			currentDBFile.readObjects();
			currentDBFile.close();
			return true;
		}
		catch (Exception e)
		{
			StsException.outputException("StsDBFile.readDatabaseObjectFile() failed.", e, StsException.WARNING);
			return false;
		}
	}

    static public boolean readObjectFile(String pathname, Object object, StsProgressBar processView)
	{
		try
		{
			StsFile file = StsFile.constructor(pathname);
			StsDBFileObject currentDBFile = StsDBFileObject.openRead(file, processView);
            if(currentDBFile == null) return false;
            StsSimpleObjectAddCmd.setObject(object);
            boolean readOk = currentDBFile.readObjects();
			currentDBFile.close();
			return readOk;
		}
		catch (Exception e)
		{
			StsException.outputException("StsDBFile.readObjectFile() failed.", e, StsException.WARNING);
			return false;
		}
	}
	protected String getDBTypeName()
	{
		return "S2S-DB-OBJECT";
	}

	//
	// TEST CODE
	//

	static private void testCursor3d(String pathname)
	{
		Main.isDbDebug = true;
		// create a set of test properties and write them out
		StsGLPanel3d glPanel3d = new StsGLPanel3d();
		StsCursor3d cursor3d = new StsCursor3d();
//		writeObjectFile(pathname, cursor3d);
	}

	static private void testProperties(String pathname)
	{
		Main.isDbDebug = true;
		// create a set of test properties and write them out
		StsProperties properties = new StsProperties();
		properties.set("integer", 1);
		properties.set("boolean", true);
		properties.set("string", "STRING");
		writeObjectFile(pathname, properties, null);

		// create an empty set of properties and populate them from the input file
		properties = new StsProperties();
		readObjectFile(pathname, properties, null);

		// check that we have the right answers
		int i = properties.getInteger("integer");
		boolean b = properties.getBoolean("boolean");
		String s = properties.getProperty("string");
		if (debug)
			System.out.println("integer: " + i + " boolean: " + b + " string: " + s);
	}

	static private void testPersistManager(String pathname)
	{
		Main.isDbDebug = true;
		StsWin3d win3d = new StsWin3d();
		StsModel model = new StsModel();
		StsSeismicVolume volume = new StsSeismicVolume(true);
		win3d.addSelectedObject(volume);
		writeObjectFile(pathname, win3d, null);
		win3d = new StsWin3d();
		readObjectFile(pathname, win3d, null);
		volume = (StsSeismicVolume)win3d.getSelectedObject(StsSeismicVolume.class);
		if (debug)
			System.out.println("Retrieved line2d[0] with index " + volume.getIndex());
	}

	static private void testCustomSerializable(String pathname)
	{
		Main.isDbDebug = true;
		StsPoint point = new StsPoint(1, 2, 3);
		StsGridSectionPoint gridPoint = new StsGridSectionPoint(point, 0, 0, true);
		writeObjectFile(pathname, gridPoint, null);
		gridPoint = new StsGridSectionPoint();
		readObjectFile(pathname, gridPoint, null);
	}

	static private void testProject(String pathname)
	{
		Main.isDbDebug = true;
		StsModel model = new StsModel();
		model.setProject(new StsProject());
		StsSeismicVolumeClass seismicVolumeClass = new StsSeismicVolumeClass();
		writeObjectFile(pathname, model, null);
		model = new StsModel();
		readObjectFile(pathname, model, null);
	}

	static private void testSegyFormat(String pathname)
	{
        StsSEGYFormat segyFormat = new StsSEGYFormat();
        segyFormat.setName("segyFormatTest");
		writeObjectFile(pathname, segyFormat, null);
        segyFormat = new StsSEGYFormat();
        readObjectFile(pathname, segyFormat, null);
        System.out.println(segyFormat.getName());
    }

    static private void testSemblanceProperties(String pathname)
	{
		Main.isDbDebug = true;
		StsModel model = new StsModel();
		model.setProject(new StsProject());
		StsSeismicVolumeClass seismicVolumeClass = new StsSeismicVolumeClass();
	}

	static private void testInnerClass(String pathname)
	{
		InnerClassTest test = new InnerClassTest();
		if (debug)
			System.out.println("in    a: " + test.a + " b: " + test.inner.b);
		writeObjectFile(pathname, test, null);
		readObjectFile(pathname, test, null);
		if (debug)
			System.out.println("out   a: " + test.a + " b: " + test.inner.b);
	}

	static public void main(String[] args)
	{
		// setup the model for testing
		StsModel model = StsModel.constructor();
		// create an object to copy
		StsObjectTest object = new StsObjectTest();
		object.initialize();
		object.setName("original");
		// output original before copy to compare with after copy
		object.print();
		String pathname = "c:\\stsdev\\c80_54\\dbtest";
		StsDBFileObject.writeObjectFile(pathname, object, null);
		StsModel newModel = StsModel.constructor();
        object = new StsObjectTest();
        StsDBFileObject.readObjectFile(pathname, object, null);
        StsObject[] stsObjects = newModel.getObjectList(StsObjectTest.class);
        for(StsObject stsObject : stsObjects)
            ((StsObjectTest)stsObject).print();
        
//        testSegyFormat("c:\\stsdev\\c80_17-tom\\segyFormat.test");
//		testProject("c:\\stsdev\\c75k-tom\\project.test");
//		testCustomSerializable("c:\\properties.test");
//		testCustomSerializable("c:\\stsdev\\c75k-tom\\properties.test");
//		dbFieldTest("c:\\dbFieldTest");
//		dbFieldTest("c:\\stsdev\\c75k-tom\\dbFieldTest2");
//		testInnerClass("c:\\stsdev\\c63c-jogl\\innertest");  // Inner classes don't work for STS persistence
//		testCursor3d("c:\\stsdev\\c63c-jogl-persist\\cursor3dtest");
//		testSemblanceProperties("c:\\stsdev\\c75k-tom\\semblanceProperties.test");
//		testProperties("c:\\stsdev\\c63c-jogl-persist\\properties.test");
//		testPersistGLPanel("c:\\stsdev\\c75k-tom\\glPanel.test");
//		testPersistManager("c:\\stsdev\\c75k-tom\\windows.test");
	}

	static private void dbFieldTest(String pathname)
	{
		Main.isDbDebug = true;
		DBFieldTest test = new DBFieldTest();
		test.initialize();
		test.print("in: ");
		writeObjectFile(pathname, test, null);
		test = new DBFieldTest();
		readObjectFile(pathname, test, null);
		test.print("out: ");
	}
}

// inner classes just don't serialize easily as they have a field reference "this$0"
// which is the outer class and they don't have a direct default constructor
class InnerClassTest extends StsSerialize
{
	int a = 1;
	InnerClass inner;
	InnerClassTest()
	{
		inner = new InnerClass();
	}

	class InnerClass implements StsSerializable
	{
		int b = 2;

		InnerClass()
		{
		}
	}
}

class DBFieldTest extends StsSerialize
{
	java.awt.Color red;
	int[][] int2null;
	int[][][] int3;
	Object[] args;
	int i = 0;

	DBFieldTest()
	{
	}

	void initialize()
	{
		int3 = new int[][][] { { { 1, 2 }, { 11, 12, 13 } },
							   { { 21, 22, 23, 24 }, { 31, 32 }, {1, 2, 3} } };
		args = new Object[] { int3 };
		i = 99;
		red = new java.awt.Color(100, 50, 25);
	}

	void print(String string)
	{
		System.out.println(string + " red = " + red.getRed());
		System.out.println(string + " green = " + red.getGreen());
		System.out.println(string + " blue = " + red.getBlue());
		printArray(string + " int2null", int2null);
		printArray(string + " int3", int3);
		int[][][] argInt3 = (int[][][])args[0];
		printArray(string + " argInt3", argInt3);
	}

	void printArray(String string, Object arrayObject)
	{
		if (arrayObject == null)
			System.out.println(string + " null");
		else
		{
			System.out.println(string);
			printArrayValues(arrayObject);
		}
	}

	void printArrayValues(Object arrayObject)
	{
//		System.out.println(string);
		if (arrayObject.getClass().isArray())
		{
			int length = Array.getLength(arrayObject);
			for (int n = 0; n < length; n++)
				printArrayValues(Array.get(arrayObject, n));
			System.out.println();
		}
		else
			System.out.print(arrayObject + "  ");
	}
}
