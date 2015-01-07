package com.Sts.DB;

import com.Sts.DB.DBCommand.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

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

public class StsDBFileObjectTrader extends StsDBFile
{
	private Collection uniqueStsObjects = new TreeSet();
	private StsObject shallowCopyObject = null;

	// *****************************************************************************************************************
	// Static Public Methods -- these represent the main API for this class
	// *****************************************************************************************************************

	static public boolean exportStsObject(String pathname, StsObject object)
	{
		return exportStsObject(pathname, object, null);
	}

	static public boolean exportStsObject(String pathname, StsObject object, StsProgressBar processView)
	{
		try
		{
			StsFile file = StsFile.constructor(pathname);
			StsDBFileObjectTrader objectTrader = StsDBFileObjectTrader.openWrite(file, processView);
			boolean result = objectTrader.deepExport(object);
			objectTrader.close();
			return result;
		}

		catch (Exception e)
		{
			new StsMessage(null, StsMessage.ERROR, "StsObjectTrader.export() failed. File: " + pathname);
			return false;
		}
	}

	static public StsObject importStsObject(String pathname)
	{
		return importStsObject(pathname, null);
	}

	static public StsObject importStsObject(String pathname, StsProgressBar processView)
	{
		try
		{
			StsFile file = StsFile.constructor(pathname);
			StsDBFileObjectTrader objectTrader = StsDBFileObjectTrader.openRead(file, processView);
			objectTrader.readObjects();
			objectTrader.close();
			return StsExportRootObjectCmd.getRootObject();
		}
		catch (Exception e)
		{
			new StsMessage(null, StsMessage.ERROR, "StsObjectTrader.export() failed. File: " + pathname);
			return null;
		}
	}

	static public StsObject clone(StsObject object, StsProgressBar processView)
	{
		try
		{
            return (StsObject)object.clone();
		}
		catch (Exception e)
		{
			new StsMessage(null, StsMessage.ERROR, "StsObjectTrader.clone() failed.");
			return null;
		}
	}

	// *****************************************************************************************************************
	// Protected Methods -- these should never need to be called by an applications developer
	// *****************************************************************************************************************

	protected String getDBTypeName()
	{
		return "S2S-DB-OBJECT-TRADER";
	}

	protected boolean add(StsObject obj)
	{
		// wrapping object to avoid class cast exceptions on comparing StsObjects directly
		return uniqueStsObjects.add(new StsObjectWrapper(obj));
	}

	protected StsDBInputStream getDBInputStream(InputStream in, StsDBInputStream oldDBInputStream) throws IOException
	{
		if (shallowCopyObject != null)
		{
			return new StsShallowImportDBInputStream(this, in, oldDBInputStream, shallowCopyObject);
		}
		else
		{
			return new StsDeepImportDBInputStream(this, in, oldDBInputStream);
		}
	}

	// *****************************************************************************************************************
	// Private Methods
	// *****************************************************************************************************************

	private StsDBFileObjectTrader(StsAbstractFile file, StsProgressBar progressBar) throws StsException, IOException
	{
		super(file, progressBar);
	}

	static private StsDBFileObjectTrader openWrite(StsAbstractFile file, StsProgressBar progressBar)
	{
		try
		{
			StsDBFileObjectTrader dbFile = new StsDBFileObjectTrader(file, progressBar);
			if (!dbFile.openWrite())
			   return null;
			return dbFile;
		}

		catch (Exception e)
		{
			new StsMessage(null, StsMessage.ERROR, "StsDBFileObjectTrader.openWrite() failed. File: " + file.getFilename());
			return null;
		}
	}

	static private StsDBFileObjectTrader openRead(StsAbstractFile file, StsProgressBar progressBar)
	{
		try
		{
			StsDBFileObjectTrader dbFile = new StsDBFileObjectTrader(file, progressBar);
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

	private boolean deepExport(StsObject object) throws IllegalAccessException
	{
		StsDBTypeObject dbType = (StsDBTypeObject)getCurrentDBType(object);
		dbType.exportObject(this, object);
		StsExportRootObjectCmd rootObjectCmd = new StsExportRootObjectCmd(object);

		// now write out all the export commands
		ArrayList commands = new ArrayList();
		Iterator itr = uniqueStsObjects.iterator();
		while (itr.hasNext())
		{
			StsObjectWrapper obj = (StsObjectWrapper)itr.next();
			if (debug) System.out.println("Writing export command for " + obj.object);
			StsInstanceExportCmd command = new StsInstanceExportCmd(obj.object);
			commands.add(command);
		}
		commands.add(rootObjectCmd);
		return commitCmdList("Export", commands);
	}

	private boolean shallowExport(StsObject object) throws IllegalAccessException
	{
		StsExportRootObjectCmd rootObjectCmd = new StsExportRootObjectCmd(object);
		StsInstanceExportCmd command = new StsInstanceExportCmd(object);
		ArrayList commands = new ArrayList();
		commands.add(command);
		commands.add(rootObjectCmd);
		return commitCmdList("Export", commands);
	}

	private class StsObjectWrapper implements Comparable
	{
		StsObject object;
		StsObjectWrapper(StsObject object)
		{
		   this.object = object;
		}

		public int compareTo(Object o)
		{
			StsObjectWrapper that = (StsObjectWrapper)o;
			return this.object.hashCode() - that.object.hashCode();
		}
	}

	class StsDeepImportDBInputStream extends StsDBInputStream
	{
		// This class remaps the indicies of all incoming objects to new indicies

		private HashMap objectsByClass = new HashMap();

		public StsDeepImportDBInputStream(StsDBFile dbFile, InputStream in, StsDBInputStream oldDBInputStream) throws IOException
		{
			super(dbFile, in, oldDBInputStream);
		}

		public StsObject getModelObject(Class c, int index)
		{
			StsObject object = getImportObject(c, index);
			if (object != null) return object;

			StsModel model = StsObject.getCurrentModel();

			StsDBTypeClass dbType = (StsDBTypeClass)getCurrentDBType(c);
			object = (StsObject)dbType.newInstance();
			model.add(object);

			storeImportObject(c, index, object);

			return object;
		}

		public StsObject getModelObjectOrNull(Class c, int index)
		{
			return getImportObject(c, index);
		}

		private StsObject getImportObject(Class c, int index)
		{
			HashMap objectsByIndex = (HashMap)objectsByClass.get(c);
			if (objectsByIndex == null)
			{
				return null;
			}
			return (StsObject)objectsByIndex.get(new Integer(index));
		}

		private void storeImportObject(Class c, int index, StsObject object)
		{
			HashMap objectsByIndex = (HashMap)objectsByClass.get(c);
			if (objectsByIndex == null)
			{
				objectsByIndex = new HashMap();
				objectsByClass.put(c, objectsByIndex);
			}
			objectsByIndex.put(new Integer(index), object);
		}
	}

	class StsShallowImportDBInputStream extends StsDBInputStream
	{
		// This class remaps only the object being copied. All other objects are read by reference to their original indicies.

		private Class mappedObjectClass = null;
		private int mappedObjectIndex = -1;
		private StsObject mappedObject = null;

		public StsShallowImportDBInputStream(StsDBFile dbFile, InputStream in, StsDBInputStream oldDBInputStream, StsObject shallowCopyObject) throws IOException
		{
			super(dbFile, in, oldDBInputStream);
			mappedObjectClass = shallowCopyObject.getClass();
			mappedObjectIndex = shallowCopyObject.getIndex();
		}

		public StsObject getModelObject(Class c, int index)
		{
			if (c == mappedObjectClass && index == mappedObjectIndex)
			{
				if (mappedObject == null)
				{
					StsModel model = StsObject.getCurrentModel();
					StsDBTypeClass dbType = (StsDBTypeClass)getCurrentDBType(c);
					mappedObject = (StsObject)dbType.newInstance();
					model.add(mappedObject);
				}
				return mappedObject;
			}
			return super.getModelObject(c, index);
		}

		public StsObject getModelObjectOrNull(Class c, int index)
		{
			if (c == mappedObjectClass && index == mappedObjectIndex)
			{
				return mappedObject;
			}
			return super.getModelObjectOrNull(c, index);
		}
	}

	public static void main(String[] args)
	{
		try
		{
            StsModel model = StsModel.constructor();
			StsObject.setCurrentModel(model);
        /*
            StsObjectRefList profileList = StsObjectRefList.constructor(10, 10, "velocityProfiles", null);
            for(int n = 0; n < 5; n++)
            {
                StsVelocityProfile profile = new StsVelocityProfile(false);
                profile.addProfilePoint(1.0f, 1.0f);
                profile.addToModel();
                profileList.add(profile);
            }

            System.out.println("original profile list");
            StsVelocityProfile profile = (StsVelocityProfile)profileList.getElement(0);
            System.out.println("    velocityProfiles count = " + model.getObjectList(StsVelocityProfile.class).length);
            System.out.println("        profile point 0: " + profile.getProfilePoint(0).toString());

            StsObjectRefList deepProfileList = (StsObjectRefList)StsDBFileObjectTrader.deepCopyStsObject(profileList, null);

            System.out.println("deep copy profile list");
            StsVelocityProfile deepProfile = (StsVelocityProfile)deepProfileList.getElement(0);
            System.out.println("    velocityProfiles count = " + model.getObjectList(StsVelocityProfile.class).length);
            System.out.println("        profile point 0: " + deepProfile.getProfilePoint(0).toString());
        */
        /*
            StsVelocityProfile profile = new StsVelocityProfile(false);
            profile.addProfilePoint(1.0f, 1.0f);
            profile.addToModel();

            System.out.println("velocityProfiles count = " + model.getObjectList(StsVelocityProfile.class).length);
            System.out.println("profile point 0: " + profile.getProfilePoint(0).toString());

            System.out.println("shallow copy of profile");
            StsVelocityProfile shallowProfile = (StsVelocityProfile)StsDBFileObjectTrader.shallowCopyStsObject(profile, null);

            System.out.println("velocityProfiles count = " + model.getObjectList(StsVelocityProfile.class).length);
            System.out.println("shallow profile point 0: " + shallowProfile.getProfilePoint(0).toString());

            System.out.println("deep copy of profile");
            StsVelocityProfile deepProfile = (StsVelocityProfile)StsDBFileObjectTrader.deepCopyStsObject(profile, null);

            System.out.println("velocityProfiles count = " + model.getObjectList(StsVelocityProfile.class).length);
            System.out.println("deep profile point 0: " + deepProfile.getProfilePoint(0).toString());
        */

            StsPreStackLineSet3d lineSet3d = new StsPreStackLineSet3d();
            StsPreStackVelocityModel3d velocityModel = new StsPreStackVelocityModel3d(lineSet3d);
            StsObjectRefList profileList = StsObjectRefList.constructor(10, 10, "velocityProfiles", velocityModel);
            for(int n = 0; n < 5; n++)
            {
                StsVelocityProfile profile = new StsVelocityProfile(false);
                profile.addProfilePoint(1.0f, 1.0f);
                profile.addToModel();
                profileList.add(profile);
            }

            System.out.println("lineSet count = " + model.getObjectList(StsPreStackLineSet3d.class).length);
            System.out.println("velocityModel count = " + model.getObjectList(StsPreStackVelocityModel3d.class).length);
            System.out.println("objRefList count = " + model.getObjectList(StsObjectRefList.class).length);
            System.out.println("velocityProfiles count = " + model.getObjectList(StsVelocityProfile.class).length);

            debugPrintVelocityModel("originalVelocityModel", velocityModel);

            StsPreStackLineSet3d shallowCopyLineSet3d = (StsPreStackLineSet3d)StsDBFileObjectTrader.clone(lineSet3d, null);

            // make a shallow copy of the velocity model and print debug
            StsPreStackVelocityModel3d shallowCopyVelocityModel = (StsPreStackVelocityModel3d)StsDBFileObjectTrader.clone(velocityModel, null);
            shallowCopyLineSet3d.setVelocityModel(shallowCopyVelocityModel);

            System.out.println("lineSet count = " + model.getObjectList(StsPreStackLineSet3d.class).length);
            System.out.println("velocityModel count = " + model.getObjectList(StsPreStackVelocityModel3d.class).length);
            System.out.println("objRefList count = " + model.getObjectList(StsObjectRefList.class).length);
            System.out.println("velocityProfiles count = " + model.getObjectList(StsVelocityProfile.class).length);

            debugPrintVelocityModel("shallowCopyVelocityModel", shallowCopyVelocityModel);

            // make a deep copy of the velocity model and print debug
            StsPreStackVelocityModel3d deepCopyVelocityModel = (StsPreStackVelocityModel3d)StsDBFileObjectTrader.clone(velocityModel, null);
            StsObjectRefList velocityProfiles = velocityModel.getVelocityProfiles();
            StsObjectRefList velocityProfilesCopy = velocityProfiles.deepCopy();
            deepCopyVelocityModel.setVelocityProfiles(velocityProfilesCopy);
            shallowCopyLineSet3d.setVelocityModel(deepCopyVelocityModel);

            System.out.println("lineSet count = " + model.getObjectList(StsPreStackLineSet3d.class).length);
            System.out.println("velocityModel count = " + model.getObjectList(StsPreStackVelocityModel3d.class).length);
            System.out.println("objRefList count = " + model.getObjectList(StsObjectRefList.class).length);
            System.out.println("velocityProfiles count = " + model.getObjectList(StsVelocityProfile.class).length);

            debugPrintVelocityModel("deepCopyVelocityModel", deepCopyVelocityModel);

        /*
            StsWellPlan newWellPlan;
            
            StsWellPlan wellPlan = new StsWellPlan("Well plan", true);
			StsPlatform platform = new StsPlatform();
			wellPlan.setPlatform(platform);

			System.out.println("Initialised model with one well plan and one platform");
			System.out.println("Well plan count = " + model.getObjectList(StsWellPlan.class).length);
			System.out.println("Platform count = " + model.getObjectList(StsPlatform.class).length);
			System.out.println("");

            System.out.println("export/import from file");
            // Export to a file
			StsDBFileObjectTrader.exportStsObject("c:\\stsdev\\c75k-tom\\exportWellPlan", wellPlan, null);
            // Import from file
			newWellPlan = (StsWellPlan)StsDBFileObjectTrader.importStsObject("c:\\stsdev\\c75k-tom\\exportWellPlan", null);
            System.out.println("Well plan count = " + model.getObjectList(StsWellPlan.class).length);
			System.out.println("Platform count = " + model.getObjectList(StsPlatform.class).length);
			System.out.println("");

			// Copy in memory using export/import mechanism
            System.out.println("Deep copy well plan");
            newWellPlan = (StsWellPlan)StsDBFileObjectTrader.deepCopyStsObject(wellPlan, null);
            System.out.println("Well plan count = " + model.getObjectList(StsWellPlan.class).length);
			System.out.println("Platform count = " + model.getObjectList(StsPlatform.class).length);
			System.out.println("");

			// Copy in memory using export/import mechanism
            newWellPlan = (StsWellPlan)StsDBFileObjectTrader.shallowCopyStsObject(wellPlan, null);
            System.out.println("Shallow copy well plan");
            System.out.println("Well plan count = " + model.getObjectList(StsWellPlan.class).length);
			System.out.println("Platform count = " + model.getObjectList(StsPlatform.class).length);
			System.out.println("");
	    */
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

    static void debugPrintVelocityModel(String message, StsPreStackVelocityModel velocityModel)
    {
        System.out.println(message);
        System.out.println("    object index: " + velocityModel.getIndex());
        StsObjectRefList velocityProfiles = velocityModel.getVelocityProfiles();
        System.out.println("    profiles list index: " + velocityProfiles.getIndex());
        int nProfiles = velocityProfiles.getSize();
        for(int n = 0; n < nProfiles; n++)
        {
            StsVelocityProfile profile = (StsVelocityProfile)velocityProfiles.getElement(n);
            System.out.println("        profile index: " + profile.getIndex());
            StsPoint semblancePoint = profile.getProfilePoint(0);
            semblancePoint.print();
        }
    }
}
