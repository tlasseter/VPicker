package com.Sts.MVC;

import com.Sts.DB.DBCommand.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */

/** This is the central class for the management of views and view properties persistence via custom serialization
 *
 */
public class StsPropertiesPersistManager implements StsSerializable //, Runnable
{
	private StsClass[] classes;
    private StsProperties modelProperties;
	transient private StsModel model;

    static final boolean debug = false;
	static final long serialVersionUID = 1l;

	public StsPropertiesPersistManager()
	{

	}

	public StsPropertiesPersistManager(StsModel m)
	{
		setModel(m);
	}

	public void setModel(StsModel m)
	{
		model = m;
	}

	/* from a save-as */
	public void setModelAs(StsModel m)
	{
		setModel(m);
	}

	public void restore()
	{
		try
		{
            model.setClasses(classes);
            if(modelProperties != null) model.properties = modelProperties;
		}
		catch (Exception e)
		{
			new StsMessage(model.win3d, StsMessage.WARNING, "Properties restoration failed.\nError: " + e.getMessage());
		}
	}

	public void save(StsDBFile db)
	{
		try
		{
            classes = model.classList;
            modelProperties = model.properties;
            if(debug) System.out.println("Writing properties to DB.");
            // db.debugCheckWritePosition("before properties persist write");
            if(!db.commitCmd("save properties", new StsAddTransientModelObjectCmd(this, "propertiesPersistManager")))
            {
                StsException.systemError(this, "save", "failed to commit cmd to db. status: " + db.statusStrings[db.status] + " transaction:" + db.transactionTypeStrings[db.transactionType]);
                return;
            }
            StsMessageFiles.logMessage("Properties saved in db.");
        }
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}

}