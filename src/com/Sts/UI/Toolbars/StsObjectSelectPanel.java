package com.Sts.UI.Toolbars;

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.WorkflowPlugIn.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.event.*;

/**
 * Created by IntelliJ IDEA.
* User: Tom Lasseter
* Date: Jan 22, 2008
* Time: 7:18:54 PM
*/

/** This class is a JPanel with a comboBox containing a list of instances belong to this class and optionally a list of subclasses
 *  selectedObject is persistent and is recovered when the objectSelectPanel is built as a component of an StsComboBoxToolbar.
 */
public class StsObjectSelectPanel extends StsJPanel implements PopupMenuListener, StsSerializable
{
    public StsObject selectedObject = null;
    transient StsModel model;
    transient Class parentClass = null;
    transient StsWorkflowPlugIn.ComboBoxDescriptor descriptor;
    transient Class[] childClasses = null;
    transient StsToggleButton toggleButton;
    transient StsJComboBox comboBox;
    transient StsComboBoxToolbar stsComboBoxToolbar;

    public StsObjectSelectPanel()
    {
    }
    
    public StsObjectSelectPanel(StsComboBoxToolbar stsComboBoxToolbar, StsWorkflowPlugIn.ComboBoxDescriptor descriptor, StsModel model)
    {
        initialize(stsComboBoxToolbar, descriptor, model);
    }
    
    public void initialize(StsComboBoxToolbar stsComboBoxToolbar, StsWorkflowPlugIn.ComboBoxDescriptor descriptor, StsModel model)
    {
        this.stsComboBoxToolbar = stsComboBoxToolbar;
        this.descriptor = descriptor;
        this.model = model;
        String[] classNames = descriptor.classNames;
        String parentClassName = descriptor.parentClassName;
        String iconName = descriptor.selectedIconName;

        if(parentClassName == null)
        {
            if(classNames == null) return;
            parentClassName = classNames[0];
        }
        parentClass = StsToolkit.getClassForName(parentClassName);
        if(parentClass == null)
        {
            StsException.systemError("StsComboBoxToolbar failed. Couldn't load class: " + parentClassName);
            return;
        }
        StsClass stsClass = stsComboBoxToolbar.model.getCreateStsClass(parentClass);

        if(classNames != null)
        {
            childClasses = new Class[classNames.length];
            int nChildClasses = 0;
            for(int n = 0; n < classNames.length; n++)
            {
                Class childClass = StsToolkit.getClassForName(classNames[n]);
                if(childClass != null)
                    childClasses[nChildClasses++] = childClass;
            }
            if(nChildClasses < classNames.length)
                childClasses = (Class[])StsMath.trimArray(childClasses, nChildClasses);
        }

        toggleButton = new StsToggleButton(iconName, "toggle visibility of " + iconName, this, "toggleOn", "toggleOff");
        toggleButton.addIcons(descriptor.selectedIconName, descriptor.deselectedIconName);
        toggleButton.setSelected(true);
        addToRow(toggleButton);

        comboBox = new StsJComboBox();
        comboBox.setLightWeightPopupEnabled(false);
        addComboBoxListeners();
        addClassListeners();
        createComboBoxList();
        addEndRow(comboBox);
    }

    public void setVisible(boolean visible)
    {
        if(isVisible() == visible) return;
        super.setVisible(visible);
        comboBox.setVisible(visible);
    }

    private void addComboBoxListeners()
    {
        comboBox.addPopupMenuListener(this);
        comboBox.addItemListener(new ItemListener()
        {
            public void itemStateChanged(ItemEvent e)
            {
                int state = e.getStateChange();
                Object object = e.getItem();
                if(state == ItemEvent.SELECTED)
                {
                    comboBoxSelectObject(object);
                }
            }
        });
    }

    private void objectPanelComboBoxSelectObject(Object object)
    {
        comboBoxSelectObject(object, true);
    }

    private void comboBoxSelectObject(Object object)
    {
        comboBoxSelectObject(object, false);
    }

    private void comboBoxSelectObject(Object object, boolean isObjectPanelSelect)
    {
        if(selectedObject == object) return;
//        stsComboBoxToolbar.win3d.glPanel3d.setObject(object);   // Don't comment this line out -- required for multiple windows.
        selectedObject = (StsObject)object;
        objectPanelSelectObject(object);
        StsWin3dBase win3d = stsComboBoxToolbar.win3d;
        if(win3d instanceof StsWin3d || win3d instanceof StsWin3dFull)
        {
            win3d.windowFamilyViewObjectChangedAndRepaint(this, object);
            // this call may not be redundant as objectPanelSelectObject eventually calls it as well
            // if setCurrentObject checks that object is already current, no problem;
            // otherwise we might have an extra redraw
            model.setCurrentObject((StsObject)object);
            toggleButton.setSelected(true);
        }
        else
        {
            StsView currentView = win3d.getGlPanel3d().getView();
            if(currentView.viewObjectChanged(this, object))
                currentView.viewObjectRepaint(this, object);
        }
//        stsComboBoxToolbar.win3d.glPanel3d.repaint();
    }

    private void objectPanelSelectObject(Object object)
    {
        if(object == null) return;

        if(stsComboBoxToolbar.win3d.isMainWindow())
        {
            StsObjectTree objectTreePanel = stsComboBoxToolbar.model.win3d.objectTreePanel;
            if(objectTreePanel != null)
                objectTreePanel.selected((StsObject)object);
        }
    }

    public void setComboBoxItem(Object object)
    {
        comboBoxSelectObject(object);
        comboBox.setSelectedItem(object);
        toggleButton.setSelected(true);
    }

    public void objectPanelSetComboBoxItem(Object object)
    {
        objectPanelComboBoxSelectObject(object);
        comboBox.setSelectedItemNoActionEvent(object);
        toggleButton.setSelected(true);
    }

    public void toggleOn()
    {
        if(selectedObject != null) model.toggleOnCursor3dObject(selectedObject);
    }

    public void toggleOff()
    {
        if(selectedObject != null) model.toggleOffCursor3dObject(selectedObject);
    }

    public StsObject deleteComboBoxObject(StsObject object)
    {
        StsObject[] objects = getComboBoxObjects();
        int index = StsMath.arrayGetIndex(objects, object);
        if(index == -1) return null;
        int length = objects.length;
        if(length == 1)
        {
            selectedObject = null;
            setVisible(false);
            return null;
        }
        if(index > 0)
            selectedObject = objects[index-1];
        else // index == 0;
            selectedObject = objects[index+1];
        return selectedObject;
    }

    /** If instances are added or deleted to this class, they may change the visibility of this objectSelectPanel on the toolbar. */
    private void addClassListeners()
    {
        try
        {
            for(int i = 0; i < childClasses.length; i++)
            {
                StsClass stsClass = stsComboBoxToolbar.model.getCreateStsClass(childClasses[i]);
                if(stsClass == null) continue;
                addClassListeners(stsClass);
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "createComboBoxList", e);
        }
    }

    private void createComboBoxList()
    {
        StsObject[] objects = getComboBoxObjects();
        createComboBoxItems(objects);
    }

    private StsObject[] getComboBoxObjects()
    {
        StsObject[] objects = new StsObject[0];
        try
        {
            for(int i = 0; i < childClasses.length; i++)
            {
                StsClass stsClass = stsComboBoxToolbar.model.getStsClass(childClasses[i]);
                if(stsClass == null) continue;
                objects = (StsObject[])StsMath.arrayAddArray(objects, stsClass.getObjectList());
            }
            return objects;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getComboBoxObjects", e);
            return objects;
        }
    }

    private void addClassListeners(StsClass stsClass)
    {
        stsClass.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                String command = e.getActionCommand();
                StsObject object = (StsObject)e.getSource();
                if(command.equals("add"))
                {
                    comboBoxAddObject(object);
                }
                else if(command.equals("delete"))
                {
                    comboBoxDeleteObject(object);
                }
                else if(command.equals("selected"))
                {
                    stsComboBoxToolbar.comboBoxSetItem(object);
                }
                else if(command.equals("remove"))
                {
                    comboBoxDeleteObject(object);
                }
            }
        });
    }

    private void comboBoxAddObject(StsObject object)
    {
        if(object == null) return;
        if(comboBox == null) return;
        selectedObject = object;
        Runnable addObjectRunnable = new Runnable()
        {
            public void run()
            {
                DefaultComboBoxModel comboBoxModel = (DefaultComboBoxModel)comboBox.getModel();
                comboBoxModel.addElement(selectedObject);
                comboBoxModel.setSelectedItem(selectedObject);
                StsWin3dBase window = stsComboBoxToolbar.win3d;
                StsWindowFamily family = model.getWindowFamily(window);
                family.viewObjectChanged(this, selectedObject);
            /*
                if(stsComboBoxToolbar.win3d.isMainWindow())
                {
                    comboBoxModel.setSelectedItem(selectedObject);
                    stsComboBoxToolbar.win3d.glPanel3d.setObject(selectedObject);
                }
            */
                setVisible(true);
            }
        };
        StsToolkit.runLaterOnEventThread(addObjectRunnable);
    }

    private void comboBoxDeleteObject(Object object)
    {
        if(comboBox == null) return;

        DefaultComboBoxModel comboBoxModel = (DefaultComboBoxModel)comboBox.getModel();
        int index = comboBoxModel.getIndexOf(object);
        if(index == -1) return;
        comboBoxModel.removeElement(object);
        if(comboBoxModel.getSize() == 0)
            setVisible(false);
        else
        {
            if(selectedObject == object)
            {
                comboBox.setSelectedIndex(index-1);
                object = comboBox.getSelectedItem();
                model.toggleOnCursor3dObject((StsObject)object);
            }
        }
    }

    private void createComboBoxItems(StsObject[] objects)
    {
        if(objects.length == 0)
        {
            setVisible(false);
        }
        else
        {
            DefaultComboBoxModel comboModel = new DefaultComboBoxModel(objects);
            comboBox.setModel(comboModel);
            initializeSelectedItem();
            setVisible(true);
        }
    }

    public boolean hasClass(Class objectClass)
    {
        if(objectClass == parentClass)
            return true;
        if(childClasses == null) return false;
        for(int n = 0; n < childClasses.length; n++)
            if(objectClass == childClasses[n])
                return true;
        return false;
    }

    private void initializeSelectedItem()
    {
        StsClass stsClass;
        stsClass = stsComboBoxToolbar.model.getCreateStsClass(parentClass);
        if(selectedObject != null)
        {
            comboBox.setSelectedItem(selectedObject);
            return;
        }
        // selectedObject is currently null, so find one
        selectedObject = stsClass.getCurrentObject();

        if(selectedObject != null)
        {
            comboBox.setSelectedItem(selectedObject);
            return;
        }
        if(childClasses == null) return;
        for(int n = 0; n < childClasses.length; n++)
        {
            stsClass = stsComboBoxToolbar.model.getCreateStsClass(childClasses[n]);
            selectedObject = stsClass.getCurrentObject();
            if(selectedObject != null)
            {
                comboBox.setSelectedItem(selectedObject);
                return;
            }
        }
    }

    public void popupMenuWillBecomeVisible(PopupMenuEvent e)
    {
        createComboBoxList();
        initializeSelectedItem();
    }

    public void popupMenuWillBecomeInvisible(PopupMenuEvent e) { }
    public void popupMenuCanceled(PopupMenuEvent e) { }

    public boolean forViewOnly() { return true; }
    
}