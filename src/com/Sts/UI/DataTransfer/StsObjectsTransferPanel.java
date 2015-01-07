package com.Sts.UI.DataTransfer;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Aug 28, 2008
 * Time: 4:33:59 PM
 * To change this template use File | Settings | File Templates.
 */
abstract public class StsObjectsTransferPanel extends StsGroupBox
{
    /** object listening to changes on this panel */
    protected StsObjectTransferListener listener;
    protected Object[] availableObjects;
    protected int width;
    protected int  height;
    protected StsJPanel transferPanel = new StsJPanel();
    protected StsButton addBtn = new StsButton("Add >", "Add selected file to right side.", this, "addObjects");
    protected StsButton removeBtn = new StsButton("< Remove", "Removed selected file from right side.", this, "removeObjects");
    protected StsButton addAllBtn = new StsButton("Add all >", "Add all files to right side.", this, "addAllObjects");
    protected StsButton removeAllBtn = new StsButton("< Remove all", "Remove all files from right side.", this, "removeAllObjects");
    protected StsJPanel selectObjectsPanel = new StsJPanel();
    protected DefaultListModel availableListModel = new DefaultListModel();
    protected DefaultListModel selectedListModel = new DefaultListModel();
    /** objects currently on the left-side list */
    protected JList availableList = new JList(availableListModel);
    /** objects currently selected and shown on the right-side list */
    protected JList selectedList = new JList(selectedListModel);

    public StsObjectsTransferPanel()
    {       
    }

    public void initialize(String label, StsObjectTransferListener listener, int width, int height)
    {
        this.listener = listener;
        this.width = width;
        this.height = height;
        super.initialize(label, true);
        constructTransferPanel();
    }

    abstract public Object[] getAvailableObjects();
    abstract public Object[] initializeAvailableObjects();

    protected void constructTransferPanel()
    {
        transferPanel.gbc.weightx = 0.0;
        transferPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
        Dimension buttonSize = new Dimension(80, 20);
        addBtn.setPreferredSize(buttonSize);
        addAllBtn.setPreferredSize(buttonSize);
        removeBtn.setPreferredSize(buttonSize);
        removeAllBtn.setPreferredSize(buttonSize);
        transferPanel.add(addBtn);
        transferPanel.add(addAllBtn);
        transferPanel.add(removeBtn);
        transferPanel.add(removeAllBtn);

        int buttonPanelWidth = 100;
        transferPanel.setPreferredSize(new Dimension(buttonPanelWidth, height));

        int panelWidth = (width - buttonPanelWidth) / 2;
        JScrollPane availableScrollPane = new JScrollPane();
        availableScrollPane.getViewport().setPreferredSize(new Dimension(panelWidth, height));
        availableScrollPane.getViewport().add(availableList, null);
        availableScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        availableScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);

        JScrollPane selectedScrollPane = new JScrollPane();
        selectedScrollPane.getViewport().setPreferredSize(new Dimension(panelWidth, height));
        selectedScrollPane.getViewport().add(selectedList, null);
        selectedScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
        selectedScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);

        selectObjectsPanel.gbc.fill = GridBagConstraints.BOTH;
        selectObjectsPanel.gbc.weighty = 1.0;
        selectObjectsPanel.gbc.weightx = 0.5;
        selectObjectsPanel.addToRow(availableScrollPane);
        selectObjectsPanel.gbc.fill = GridBagConstraints.NONE;
        selectObjectsPanel.gbc.weightx = 0.0;
        selectObjectsPanel.addToRow(transferPanel);
        selectObjectsPanel.gbc.fill = GridBagConstraints.BOTH;
        selectObjectsPanel.gbc.weightx = 0.5;
        selectObjectsPanel.addToRow(selectedScrollPane);

        gbc.fill = GridBagConstraints.BOTH;
        gbc.weighty = 1.0;
        addEndRow(selectObjectsPanel);

        availableList.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        availableList.addListSelectionListener
        (
            new ListSelectionListener()
    {
                public void valueChanged(ListSelectionEvent e)
                {
                    availableObjectSelected(e);
                }
            }
        );
        selectedList.addListSelectionListener
        (
            new ListSelectionListener()
            {
                public void valueChanged(ListSelectionEvent e)
                {
                    selectedObjectSelected(e);
                }
            }
        );
    }

    public void selectedObjectSelected(ListSelectionEvent e)
    {
        if(e.getValueIsAdjusting()) return;
        Object source = e.getSource();
        if (!(source instanceof JList)) return;
        Object selectedObject = selectedList.getSelectedValue();
        if (selectedObject == null) return;
        if(listener != null) listener.objectSelected(selectedObject);
    }

    public void availableObjectSelected(ListSelectionEvent e)
    {
        if(e.getValueIsAdjusting()) return;
        Object source = e.getSource();
        if (!(source instanceof JList)) return;
        Object selectedObject = availableList.getSelectedValue();
        if (selectedObject == null) return;
        if(listener != null) listener.objectSelected(selectedObject);
    }

    public boolean addAvailableObjectOk(Object object) { return true; }

    public void selectSingleVolProgrammatically(String selectedString)
    {
        selectedList.setSelectedIndex(selectedListModel.indexOf(selectedString));
    }

    public void selectSingleVolProgrammatically(int index)
    {
        selectedList.setSelectedIndex(index);
    }

    public int getSelectedCount()
    {
        return selectedListModel.size();
    }

    public int[] getSelectedIndices() { return selectedList.getSelectedIndices(); }

    public Object getSelectedObject() { return selectedList.getSelectedValue(); }

    public void setAvailableObjects(Object[] availableObjects)
    {
        this.availableObjects = availableObjects;
        availableListModel.clear();
        if(availableObjects == null) return;
        for (int n = 0; n < availableObjects.length; n++)
        {
          if(listener.addAvailableObjectOk(availableObjects[n]))
            availableListModel.addElement(availableObjects[n]);
        }
    }

    public void addObjects()
    {
        int[] selectedIndices = availableList.getSelectedIndices();
        // Must run in reverse otherwise the available list indices get screwed up on multi-select.
        for (int i = selectedIndices.length-1; i >= 0; i--)
            addObject(selectedIndices[i]);
        if (availableListModel.size() > 0) availableList.setSelectedIndex(0);
    }

    private void addObject(int index)
    {
        Object object = availableListModel.getElementAt(index);
        availableListModel.removeElement(object);
        if(selectedListModel.contains(object)) return;
        selectedListModel.addElement(object);
        listener.addObjects(new Object[] { object });
    }

    public void removeObjects()
    {
        int[] selectedIndices = selectedList.getSelectedIndices();
        for (int i = selectedIndices.length-1; i >= 0; i--)
            removeObject(selectedIndices[i]);
        if (selectedListModel.size() > 0) selectedList.setSelectedIndex(0);
    }

    private void removeObject(int index)
    {
        Object object = selectedListModel.getElementAt(index);
        selectedListModel.removeElement(object);
        availableListModel.addElement(object);
        listener.removeObjects(new Object[] { object });
    }

    public void addAllObjects()
    {
        int nObjects = availableListModel.size();
        Object[] availableObjects = availableListModel.toArray();
        availableListModel.clear();
        for (int n = 0; n < nObjects; n++)
            if(!selectedListModel.contains(availableObjects[n]))
                selectedListModel.addElement(availableObjects[n]);
        if(listener != null) listener.addObjects(availableObjects);
    }

    public void removeAllObjects()
    {
        int nObjects = selectedListModel.size();
        Object[] selectedObjects = selectedListModel.toArray();
        for (int i = 0; i < nObjects; i++)
        {
            Object selectedObject = selectedObjects[i];
            availableListModel.addElement(selectedObject);
        }
        selectedListModel.clear();
        if(listener != null) listener.removeObjects(selectedObjects);
        if (availableListModel.size() > 0) availableList.setSelectedIndex(0);
    }

    public Object[] getSelectedObjects()
    {
        return selectedListModel.toArray();
    }

    public int getNSelectedObjects()
    {
        return selectedListModel.getSize();
    }

    public int getSelectedIndex() { return selectedList.getSelectedIndex(); }

    public static void main(String[] args)
    {
        TestObjectsTransferPanelListener listener = new TestObjectsTransferPanelListener();
        try
        {
            UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.WindowsLookAndFeel());
        }
        catch (Exception e) { }
        StsModel.constructor();
        StsWell[] wells = new StsWell[5];
        for(int n = 0; n < 5; n++)
            wells[n] = new StsWell("Well-" + n, true);
        TestObjectsTransferPanel panel = new TestObjectsTransferPanel(wells, listener, 400, 100);
        com.Sts.Utilities.StsToolkit.createDialog(panel);
    }
}

class TestObjectsTransferPanel extends StsObjectsTransferPanel
{
    Object[] objects;

    TestObjectsTransferPanel(Object[] objects, StsObjectTransferListener listener, int width, int height)
    {
        super.initialize("Test Objects Transfer Panel", listener, width, height);
        setAvailableObjects(objects);
    }

    public Object[] getAvailableObjects()
    {
        return objects;
    }

    public Object[] initializeAvailableObjects()
    {
        objects = new String[] { "object1", "object2", "object3" };
        return objects;
    }
}

class TestObjectsTransferPanelListener implements StsObjectTransferListener
{

    TestObjectsTransferPanelListener()
    {
    }

    public void addObjects(Object[] objects)
    {
        System.out.println("add Objects:");
        printObjects(objects);
    }

    private void printObjects(Object[] objects)
    {
        for (int n = 0; n < objects.length; n++)
            System.out.println("    " + objects[n].toString());
    }

    public void removeObjects(Object[] objects)
    {
        System.out.println("remove Objects:");
        printObjects(objects);
    }

    public void objectSelected(Object selectedObject)
    {
        System.out.println("selected Object:" + selectedObject.toString());
    }

    public boolean addAvailableObjectOk(Object object) { return true; }
}