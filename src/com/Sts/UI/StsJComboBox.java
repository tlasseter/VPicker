package com.Sts.UI;

import javax.swing.*;
import java.awt.event.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 20, 2009
 * Time: 11:25:42 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsJComboBox extends JComboBox
{
    // Flag to ensure the we don't get multiple ActionEvents on item selection.
     private boolean selectingItem = false;

    public StsJComboBox()
    {
    }

    public void setSelectedItem(Object anObject)
    {
        setSelectedItem(anObject, true);
    }

    public void setSelectedItemNoActionEvent(Object anObject)
    {
        setSelectedItem(anObject, false);
    }
    
    public void setSelectedItem(Object anObject, boolean fireActionEvent)
    {
        Object oldSelection = selectedItemReminder;
        Object objectToSelect = anObject;
        if(oldSelection == null || !oldSelection.equals(anObject))
        {
            if(anObject != null && !isEditable())
            {
                // For non editable combo boxes, an invalid selection
                // will be rejected.
                boolean found = false;
                for(int i = 0; i < dataModel.getSize(); i++)
                {
                    Object element = dataModel.getElementAt(i);
                    if(anObject.equals(element))
                    {
                        found = true;
                        objectToSelect = element;
                        break;
                    }
                }
                if(!found)
                {
                    return;
                }
            }

            // Must toggle the state of this flag since this method
            // call may result in ListDataEvents being fired.
            selectingItem = true;
            dataModel.setSelectedItem(objectToSelect);
            selectingItem = false;

            if(selectedItemReminder != dataModel.getSelectedItem())
            {
                // in case a users implementation of ComboBoxModel
                // doesn't fire a ListDataEvent when the selection
                // changes.
                selectedItemChanged();
            }
        }
        if(fireActionEvent) fireActionEvent();
    }

    public void actionPerformed(ActionEvent e)
    {
        Object newItem = getEditor().getItem();
        setPopupVisible(false);
        getModel().setSelectedItem(newItem);
        String oldCommand = getActionCommand();
        setActionCommand("comboBoxEdited");
        fireActionEvent();
        setActionCommand(oldCommand);
    }
}
