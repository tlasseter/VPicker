package com.Sts.UI;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.UI.Icons.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

public class StsToggleButtonComboBox extends JComboBox // implements ActionListener
{
    private StsToggleButton currentButton = null; // currently displayed button
    private StsToggleButton selectedButton = null;  // this button is selected (may not be current one)
    private DefaultComboBoxModel model = new DefaultComboBoxModel();
    private ToggleButtonEditor toggleButtonEditor;

    public StsToggleButtonComboBox()
    {
        setModel(model);
        setRenderer(new ToggleButtonRenderer());
        toggleButtonEditor = new ToggleButtonEditor();
        setEditor(toggleButtonEditor);
        setEditable(true);
//         this.addActionListener(this);
    }

    public StsToggleButtonComboBox(ArrayList buttons)
    {
        this();
        Object button = null;
        for(int n = 0; n < buttons.size(); n++)
        {
            currentButton = (StsToggleButton)buttons.get(n);
            model.addElement(button);
        }
        if(currentButton != null) setSelectedItem(currentButton);
    }

    public void addButton(StsToggleButton button)
    {
        if(button == null) return;
        model.addElement(button);
//        currentButton = button;
//        setSelectedItem(currentButton);
    }

    public void addButton(StsButton button)
    {
        if(button == null) return;
        model.addElement(button);
//        currentButton = button;
//        setSelectedItem(currentButton);
    }

    public void addButton(String name, String iconName)
    {
        Icon icon = StsIcon.createIcon(iconName);
        StsToggleButton button = new StsToggleButton();
        button.setText(name);
        button.setIcon(icon);
        addButton(button);
    }

    public void deleteButton(StsToggleButton button)
    {
        String name = button.getName();

        int nButtons = model.getSize();
        for(int n = 0; n < nButtons; n++)
        {
            StsToggleButton listButton = (StsToggleButton)model.getElementAt(n);
            if(name.equals(listButton.getText()))
            {
                model.removeElement(listButton);
                if(nButtons == 1)
                    this.setEnabled(false);
                else
                    this.setSelectedIndex(0);
                return;
            }
        }
    }

    public void deleteButton(String name)
    {
        int nButtons = model.getSize();
        for(int n = 0; n < nButtons; n++)
        {
            StsToggleButton button = (StsToggleButton)model.getElementAt(n);
            if(name.equals(button.getText()))
            {
                model.removeElement(button);
                if(nButtons == 1)
                    this.setEnabled(false);
                else
                    this.setSelectedIndex(0);
                return;
            }
        }
    }

    // a toggle button item has been selected: toggle it selected

    public void actionPerformed(ActionEvent actionEvent)
    {
        if(!actionEvent.getActionCommand().equals("comboBoxChanged")) return;

        StsToggleButton button = (StsToggleButton)getSelectedItem();
//        button.setSelected(true);
        toggleButtonEditor.setSelectedButton(button);
    }

    class ToggleButtonRenderer extends StsToggleButton implements ListCellRenderer
    {
        public ToggleButtonRenderer()
        {
            setOpaque(true);
        }
        public Component getListCellRendererComponent(JList list, Object value, int index, boolean isSelected, boolean cellHasFocus)
        {
            return (StsToggleButton)value;
        }
    }

    class ToggleButtonEditor implements ComboBoxEditor
    {
        StsToggleButton editorToggleButton = new StsToggleButton();
//        boolean buttonPressed;
        Border downBorder = BorderFactory.createLoweredBevelBorder();
        Border upBorder = BorderFactory.createRaisedBevelBorder();

        public ToggleButtonEditor()
        {
//            buttonPressed = false;
            setEditorButtonPressed(false);
            editorToggleButton.addMouseListener(new MouseAdapter()
            {
                public void mouseReleased(MouseEvent e)
                {
                    StsToggleButton button = null;

                    int nButtons = model.getSize();
                    for(int n = 0; n < nButtons; n++)
                    {
                        button = (StsToggleButton)model.getElementAt(n);
                        if(button.getText().equals(editorToggleButton.getText())) break;
                    }
                    if(button == null) return;

                    boolean buttonPressed = button.isSelected();
                    buttonPressed = !buttonPressed;
                    setSelectedButton(button, buttonPressed);
                /*
                    if(selectedButton != null && selectedButton != button)
                    {
                        selectedButton.setSelected(false);
                        selectedButton.fireActionPerformed();
                        selectedButton = null;
                    }

                    currentButton = button;
                    if(buttonPressed) selectedButton = button;
                    else              selectedButton = null;

                    button.setSelected(buttonPressed);
                    button.fireActionPerformed();
                    setEditorButtonPressed();
                */
                }
            });
        }

        void setSelectedButton(StsToggleButton button)
        {
            boolean buttonPressed = button.isSelected();
            setSelectedButton(button, buttonPressed);
        }

        void setSelectedButton(StsToggleButton button, boolean buttonPressed)
        {
    //        super.setSelectedItem(button);
            if(selectedButton != null && selectedButton != button)
            {
                selectedButton.setSelected(false);
                selectedButton.fireActionPerformed();
                selectedButton = null;
            }

            currentButton = button;
            if(buttonPressed) selectedButton = button;
            else              selectedButton = null;

            button.setSelected(buttonPressed);
            button.fireActionPerformed();
            setEditorButtonPressed(buttonPressed);
        }

        private void setEditorButtonPressed(boolean buttonPressed)
        {
            editorToggleButton.setSelected(buttonPressed);
            if(!buttonPressed) editorToggleButton.setBorder(upBorder);
            else               editorToggleButton.setBorder(downBorder);
        }

        public Component getEditorComponent() { return editorToggleButton; }

        public Object getItem() { return editorToggleButton; }

        public void selectAll() { } // required for ComboBoxModel interface

        public void setItem(Object value)
        {
            if(value == null) return;
            StsToggleButton toggleButton = (StsToggleButton)value;
            editorToggleButton.setText(toggleButton.getText());
            editorToggleButton.setIcon(toggleButton.getIcon());
//            buttonPressed = true;
//            buttonPressed = toggleButton.isSelected();
//            setEditorButtonPressed(true);
//            editorToggleButton.setSelected(buttonPressed);
//            editorToggleButton.setBorder(upBorder);
//            buttonPressed = false;
        }

        public void addActionListener(ActionListener listener)
        {
            listenerList.add(ActionListener.class, listener);
        }
        public void removeActionListener(ActionListener listener)
        {
            listenerList.remove(ActionListener.class, listener);
        }
        public void fireActionPerformed(ActionEvent e)
        {
            Object[] listeners = listenerList.getListenerList();
            for(int i = listeners.length-2; i>= 0; i-=2)
            {
                if(listeners[i]==ActionListener.class)
                {
                    ((ActionListener)listeners[i+1]).actionPerformed(e);
                }
            }
        }
    }

    public void setSelectedButton(StsToggleButton button, boolean buttonPressed)
    {
        toggleButtonEditor.setSelectedButton(button, buttonPressed);
    }

    static public void main(String[] args)
    {
        Test test = new Test();
        StsToggleButtonComboBox comboBox = new StsToggleButtonComboBox();

        StsToggleButton seismicButton = new StsToggleButton("Seismic", "Seismic tip.", test, "doSeismic");
        seismicButton.setText("Seismic");
        seismicButton.setIcon(StsIcon.createIcon("seismic.gif"));
        comboBox.addButton(seismicButton);

        StsToggleButton surfacesButton = new StsToggleButton("Surfaces", "Surfaces tip.", test, "doSurfaces");
        surfacesButton.setText("Surfaces");
        surfacesButton.setIcon(StsIcon.createIcon("surfaces.gif"));
        comboBox.addButton(surfacesButton);

        StsToolkit.createDialog(comboBox, true);
    }
}

class Test
{
    public Test()
    {
    }

    public void doSeismic()
    {
        System.out.println("Do seismic.");
    }

    public void doSurfaces()
    {
        System.out.println("Do surfaces.");
    }
}
