//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.MVC;

import com.Sts.DB.DBCommand.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;
import com.Sts.WorkflowPlugIn.*;

import java.util.*;

/**
 *
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2006</p>
 * <p>Company: S2S Systems LLC</p>
 * @author JB West
 * @version beta 1.0
 *
 */

/**
 * This is the central class for the management of views and view properties persistence
 * Information is written to a separate database in the project binary file directory
 * with the name project-name.vws.
 */
public class StsViewPersistManager implements StsSerializable
{
    /** List of window families (parent and any number of controlled children) */
    public StsWindowFamily[] families = null;
    /** plugin * */
    String workflowName;
    /** well views */
    StsWellFrameViewModel[] wellViewModels;
    /** the principle 3d window where action takes place */
    transient public StsWin3d win3d = null;
    transient private StsModel model;
//	transient public String objectFilename = null;
//	transient private StsWindowActionManager actionManager;
//	transient private boolean initialPanelConfiguration;

    static final long serialVersionUID = 1l;
    static final boolean debug = false;

    public StsViewPersistManager()
    {
    }

    public StsViewPersistManager(StsModel m)
    {
        setModel(m);
    }

    public void setModel(StsModel m)
    {
        model = m;
//		String dirname = model.project.getRootDirString();
//		String filename = "vws." + model.getName();
//		objectFilename = new String(dirname + File.separator + "BinaryFiles" + File.separator + filename);

//		families = new ArrayList();
    }

    /* from a save-as */
    public void setModelAs(StsModel m)
    {
        model = m;
//		String dirname = model.project.getRootDirString();
//		String filename = "vws." + model.getName();
//		objectFilename = new String(dirname + File.separator + "BinaryFiles" + File.separator + filename);
    }

    /*
        public StsWin3dBase newStsWin3dBase(StsModel model, StsWin3dBase other)
        {
            return new StsWin3dBase(model, other);
        }
    */
    public void save(StsDBFile db)
    {
        if (model == null) return;
//        isGridCoordinates = model.getBooleanProperty("isGridCoordinates");
//        int nWindows = 0;
        try
        {
            for (int fam = 0; fam < families.length; fam++)
            {
                StsWindowFamily family = families[fam];
                StsWin3dBase[] windows = family.getWindows();
                for (int iwin = 0; iwin < windows.length; iwin++)
                {
                    windows[iwin].saveGeometry();
                    //                nWindows++;
                }

                StsViewTimeSeries[] timeWin = family.getTimeWindows();
                if(timeWin != null)
                {
                    for (int iwin = 0; iwin < timeWin.length; iwin++)
                    {
                        StsViewTimeSeries timeWindow = timeWin[iwin];
                        //timeWindow.                                  // Need a saveGeometry Method for time windows.
                    }
                }

                StsViewGauges[] gaugeWin = family.getGaugeWindows();
                if(gaugeWin != null)
                {
                    for (int iwin = 0; iwin < gaugeWin.length; iwin++)
                    {
                        StsViewGauges gaugeWindow = gaugeWin[iwin];
                        //                                  // Need a saveGeometry Method for gauge windows.
                    }
                }
            }
            //		saveToolbarsState(nWindows);
            saveWellViews();

            StsWorkflowPlugIn plugin = model.getWorkflowPlugIn();
            Class pluginClass = plugin.getClass();
            // remember the workflow plugin class name to reload
            workflowName = pluginClass.getName();
            // db.debugCheckWritePosition("before view persist write");
            if (!db.commitCmd("save views", new StsAddTransientModelObjectCmd(this, "viewPersistManager")))
            {
                StsException.systemError(this, "save", "failed to commit cmd to db. status: " + db.statusStrings[db.status] + " transaction:" + db.transactionTypeStrings[db.transactionType]);
                return;
            }
            // db.debugCheckWritePosition("after view persist write");
            /*
                File f = new File(objectFilename);
                if(f.exists())f.delete();
                if(debug) System.out.println("Writing views DB to file: " + objectFilename);
                StsDBFileObject.writeObjectFile(objectFilename, this, null);
    //			StsParameterFile.writeSerializedFields(this, parameterFileName);
            */
            StsMessageFiles.logMessage("Views settings saved in db");
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "save", e);
        }
    }

    private void saveWellViews()
    {
        StsObject[] wellObjects = model.getObjectList(StsWell.class);
        int length = wellObjects.length;
        if (length == 0) return;

        wellViewModels = new StsWellFrameViewModel[length];

        try
        {
            int nModels = 0;
            for (int n = 0; n < length; n++)
            {
                StsWell well = (StsWell) wellObjects[n];
                StsWellFrameViewModel wellViewModel = well.wellViewModel;
                if (wellViewModel == null) continue;
                wellViewModels[nModels++] = wellViewModel;
            }
            wellViewModels = (StsWellFrameViewModel[]) StsMath.trimArray(wellViewModels, nModels);
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "saveWellViews", e);
        }
    }

    /*
        public void run()
        {
            // not threadsafe !
            win3d = constructMainWindow(model, actionManager, initialPanelConfiguration);
    //		win3d = new StsWin3d(model, actionManager, initialPanelConfiguration);
        }
    */
    public StsWin3d newStsWin3d(StsModel model, boolean initialPanelConfiguration)
    {
        this.model = model;
        runConstructMainWindow(initialPanelConfiguration);
        return win3d;
    }

    /** read views database and construct views, restoring positions and states */
    public StsWin3d restore(boolean initialPanelConfiguration)
    {
        if (workflowName != null) model.loadWorkflowPlugIn(workflowName);

        if (families == null || families.length == 0)
            win3d = runConstructMainWindow(initialPanelConfiguration);
        else
        {
            win3d = (StsWin3d) families[0].getParent();
            initializeWindowsTransients();
        }
        restoreWellViews();
        return win3d;
    }

    /**
     * initialize transients for windows and components; defer any actions which involve window interaction until
     * subsequent startWindows() method.
     *
     * @return true if initialized ok
     */
    private boolean initializeWindowsTransients()
    {
        try
        {
            for (int fam = 0; fam < families.length; fam++)
            {
                StsWindowFamily family = families[fam];
                if (family == null) continue;
                StsWin3dBase[] win = family.getWindows();
                if (win == null) continue;
                StsWin3dFull parentWindow = family.getParent();
                if (parentWindow == null) continue;
                for (int iwin = 0; iwin < win.length; iwin++)
                {
                    StsWin3dBase window = win[iwin];
                    window.initializeTransients(model, parentWindow);
                }
                StsViewTimeSeries[] timeWin = family.getTimeWindows();
                if (timeWin != null)
                {
                    for (int iwin = 0; iwin < timeWin.length; iwin++)
                    {
                        StsViewTimeSeries timeWindow = timeWin[iwin];
                        timeWindow.initializeTransients(parentWindow, model);
                    }
                }
                StsViewGauges[] gaugeWin = family.getGaugeWindows();
                if (gaugeWin != null)
                {
                    for (int iwin = 0; iwin < gaugeWin.length; iwin++)
                    {
                        StsViewGauges gaugeWindow = gaugeWin[iwin];
                        gaugeWindow.initializeTransients(parentWindow, model);
                    }
                }
            }
            StsWin3dFull[] parentWindows = model.getParentWindows();
            for (StsWin3dFull parentWindow : parentWindows)
                parentWindow.checkAddDisplayableSections();

            Iterator<StsWin3dBase> windowIterator = model.getFamilyWindowIterator();
            while (windowIterator.hasNext())
            {
                StsWin3dBase window = windowIterator.next();
                window.start();
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "initializeWindowsTransients", e);
            return false;
        }
    }

    /*
        private boolean startWindows()
        {
             try
            {
                for(int fam = 0; fam < families.length; fam++)
                {
                    StsWindowFamily family = families[fam];
                    if(family == null) continue;
                    StsWin3dBase[] window = family.getWindows();
                    if(window == null) continue;
                    StsWin3dBase parentWindow = family.windows[0];
                    if(parentWindow == null) continue;
                    parentWindow.startWindow();
                    for(int i = 1; i < window.length; i++)
                        window[i].startWindow();

                    StsViewTimeSeries[] timeWindows = family.getTimeWindows();
                    if(timeWindows != null)
                    {
                        for (int iwin = 0; iwin < timeWindows.length; iwin++)
                            timeWindows[iwin].start();
                    }
                }
                return true;
            }
            catch(Exception e)
            {
                StsException.outputWarningException(this, "startWindows", e);
                return false;
            }
        }
    */
    public StsWin3d constructDefaultMainWindow()
    {
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                win3d = StsWin3d.constructDefaultMainWindow(model);
                model.win3d = win3d;
                StsGLPanel3d glPanel3d = win3d.getGlPanel3d();
                if (glPanel3d != null) model.setGlPanel3d(glPanel3d);
            }
        };
        StsToolkit.runWaitOnEventThread(runnable);
        return win3d;
    }

    private StsWin3d runConstructMainWindow(boolean _initialPanelConfiguration)
    {
        final boolean initialPanelConfiguration = _initialPanelConfiguration;
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                if (initialPanelConfiguration)
                    win3d = StsWin3d.constructSplashMainWindow(model);
                else
                    win3d = StsWin3d.constructDefaultMainWindow(model);

                model.win3d = win3d;
                // TODO: remove references to a "main" glPanel3d from model
                StsGLPanel3d glPanel3d = win3d.getGlPanel3d();
                if (glPanel3d != null) model.setGlPanel3d(glPanel3d);
            }
        };
        StsToolkit.runWaitOnEventThread(runnable);
        return win3d;
    }

    private void restoreWellViews()
    {
        if (wellViewModels == null) return;
        try
        {
            for (int n = 0; n < wellViewModels.length; n++)
                wellViewModels[n].restoreWellView();
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "restoreWellViews", e);
        }
    }

    StsWindowFamily newWindowFamily(StsWin3d parent)
    {
        return new StsWindowFamily(parent);
    }

    StsWindowFamily newWindowFamily(StsWin3dFull parent)
    {
        return new StsWindowFamily(parent);
    }

    public void addFamily(StsWindowFamily family)
    {
        families = (StsWindowFamily[]) StsMath.arrayAddElement(families, family);
    }

    public void deleteFamily(StsWindowFamily family)
    {
        families = (StsWindowFamily[]) StsMath.arrayDeleteElement(families, family);
    }

    public StsWin3dBase[] getWindows()
    {
        int nWindows = 0;
        for (int n = 0; n < families.length; n++)
            nWindows += families[n].windows.length;
        StsWin3dBase[] windows = new StsWin3dBase[nWindows];
        int nWindow = 0;
        for (int n = 0; n < families.length; n++)
        {
            int nFamilyWindows = families[n].windows.length;
            for (int i = 0; i < nFamilyWindows; i++)
                windows[nWindow++] = families[n].windows[i];
        }
        return windows;
    }
}