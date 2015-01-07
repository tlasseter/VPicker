package com.Sts.UI.Toolbars;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;
import com.Sts.WorkflowPlugIn.*;

import javax.swing.*;
import javax.swing.border.*;
import java.util.*;

/** ComboBoxToolbar is actually a set of comboBoxes.  Each comboBox is a set of instances of a particular group.
 *  When an instance of an stsClass belonging to this group is instantiated, the comboBox is created and this instance is added.
 *  See ObjectSelectPanel (innerclass) for the details of the comboBox.
 */
public class StsComboBoxToolbar extends StsToolbar implements StsSerializable
{
    transient public StsObjectSelectPanel[] objectSelectPanels;
    transient StsModel model = null;
    transient StsWin3dBase win3d = null;
    transient ClassLoader classLoader;

    static public final String NAME = "Object Selection Toolbar";
    static public final boolean defaultFloatable = true;

    static final String NONE = "none";
    static final String[] noneString = new String[] { NONE };

    static final boolean debug = false;

    public StsComboBoxToolbar()
    {
    }

    private StsComboBoxToolbar(StsWin3dBase win3d)
    {
        super(NAME);
        initialize(win3d);
    }

    public boolean initialize(StsWin3dBase win3d)
    {
        this.win3d = win3d;
        model = win3d.getModel();
        ArrayList comboBoxDescriptors = model.getComboBoxToolbarPlugInDescriptors();
        if(comboBoxDescriptors.size() == 0) return false;
        setName(NAME);
        setFloatable(true);

        Border border = BorderFactory.createEtchedBorder();
        setBorder(border);

        classLoader = getClass().getClassLoader();
        int nDescriptors = comboBoxDescriptors.size();

        objectSelectPanels = new StsObjectSelectPanel[nDescriptors];
        for(int n = 0; n < nDescriptors; n++)
            objectSelectPanels[n] = new StsObjectSelectPanel();

        for(int n = 0; n < nDescriptors; n++)
        {
            StsWorkflowPlugIn.ComboBoxDescriptor descriptor = (StsWorkflowPlugIn.ComboBoxDescriptor)comboBoxDescriptors.get(n);
            if(descriptor.classNames == null) continue;
            objectSelectPanels[n].initialize(this, descriptor, model);
            add(objectSelectPanels[n]);
        }
        addSeparator();
        addCloseIcon(win3d);
        return true;
    }

    /** This is a unit test initialization, not used in production; see main() method in this class */
    public boolean initialize(StsModel model, String[] classnames)
    {
        int nDescriptors = classnames.length;

        objectSelectPanels = new StsObjectSelectPanel[nDescriptors];
        for(int n = 0; n < nDescriptors; n++)
        {
            String classname = classnames[n];
            StsWorkflowPlugIn.ComboBoxDescriptor descriptor = new StsWorkflowPlugIn.ComboBoxDescriptor(classname, classname);
            objectSelectPanels[n] = new StsObjectSelectPanel(this, descriptor, model);
            add(objectSelectPanels[n]);
        }
        addSeparator();
        addCloseIcon(win3d);
        return true;
    }

    static public StsComboBoxToolbar constructor(StsWin3dBase win3d)
    {
        StsComboBoxToolbar toolbar = new StsComboBoxToolbar(win3d);
        return toolbar;
    }

    /*
    *  Called when comboBox item has been selected. Propagate this change
    *  to other interested parties and repaint.
    */

    private StsObjectSelectPanel getComboBoxPanel(Object object)
    {
        Class objectClass = object.getClass();
        for(int n = 0; n < objectSelectPanels.length; n++)
        {
            StsObjectSelectPanel objectSelectPanel = objectSelectPanels[n];
            if(objectSelectPanel != null && objectSelectPanel.hasClass(objectClass))
                return objectSelectPanel;

        }
        return null;
    }

    public StsObject getSelectedObject()
    {
        for(int n = 0; n < objectSelectPanels.length; n++)
        {
            StsObjectSelectPanel objectSelectPanel = objectSelectPanels[n];
            if(objectSelectPanel != null)
            {
                if(objectSelectPanel.selectedObject != null)
                    return objectSelectPanel.selectedObject;
            }
        }
        return null;
    }
    /**
     * Called when a comboBox item has been changed externally and we need to
     * change the item displayed in the comboBox.  By changing the model, we do
     * not fire an item stateChanged causing an endless loop.
     */
    public void comboBoxSetItem(Object object)
    {
        if(object == null) return;
        StsObjectSelectPanel objectSelectPanel = getComboBoxPanel(object);
        if(objectSelectPanel == null) return;
        objectSelectPanel.setComboBoxItem(object);
    }

    public void objectPanelSetComboBoxItem(Object object)
    {
        if(object == null) return;
        StsObjectSelectPanel objectSelectPanel = getComboBoxPanel(object);
        if(objectSelectPanel == null) return;
        objectSelectPanel.objectPanelSetComboBoxItem(object);
    }
    public StsObject deleteComboBoxObject(StsObject object)
    {
        if(object == null) return null;
        StsObjectSelectPanel objectSelectPanel = getComboBoxPanel(object);
        if(objectSelectPanel == null) return null;
        return objectSelectPanel.deleteComboBoxObject(object);

    }

    static public void main(String[] args)
    {
        try
        {
            StsModel model = StsModel.constructor();
            String objectFilename = "c:\\scratch";
            StsComboBoxToolbar toolbar = new StsComboBoxToolbar();
            String[] classnames = new String[] { StsSeismicVolume.class.getName(), StsVirtualVolume.class.getName() };
            toolbar.initialize(model, classnames);
            com.Sts.MVC.Main.isDbIODebug = true;
            StsDBFileObject.writeObjectFile(objectFilename, toolbar, null);           
            StsComboBoxToolbar newToolbar = new StsComboBoxToolbar();
            boolean success = StsDBFileObject.readDatabaseObjectFile(objectFilename, newToolbar, null);
        }
        catch(Exception e)
        {
            StsException.outputException("StsObjectSelectPanel.main() failed.", e, StsException.WARNING);
        }
    }
    
    public boolean forViewOnly()
    {
        return true;
    }

    // if added, only update the view in the main window



    /*
        private void comboBoxRemove(Object object)
        {
            JPanel comboBoxPanel = getComboBoxPanel(object);
            if (comboBoxPanel == null) return;

            JComboBox comboBox = (JComboBox) comboBoxPanel.getComponent(1);
            if (comboBox == null) return;

            DefaultComboBoxModel comboBoxModel = (DefaultComboBoxModel) comboBox.getModel();
            comboBoxModel.removeElement(object);
            comboBox.setVisible(false);

            StsClass stsClass = model.getCreateStsClass(object.getClass());
            if (stsClass == null) return;

            int size = stsClass.getSize();
            if (size == 0)
            {
                comboBox.setVisible(false);
                comboBoxPanel.setVisible(false);
                Iterator iter = objectSelectPanels.values().iterator();
                boolean visible = false;
                while (iter.hasNext())
                {
                    comboBoxPanel = (JPanel) iter.next();
                    JComboBox box = (JComboBox) comboBoxPanel.getComponent(1);
                    if (box.isVisible())
                    {
                        visible = true;
                        break;
                    }
                }
                if (!visible)
                {
                    this.setVisible(false);
                }
            }
            else
            {
                comboBox.setSelectedIndex(0);
                object = comboBox.getSelectedItem();
                win3d.glPanel3d.cursor3d.setObject(object);
            }
        }
    */

}
