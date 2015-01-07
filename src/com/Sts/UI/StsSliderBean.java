//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.UI.Beans.*;
import com.Sts.UI.Icons.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

public class StsSliderBean extends StsJPanel implements ItemListener, ChangeListener
{
    private StsJSlider slider = new StsJSlider();
    private int sliderMax = 10;
    /** if programmatically locked, this slider will not respond */
    private boolean locked = false;
//    DecimalFormat valueFormat = new DecimalFormat("###,##0.00");
//    DecimalFormat incFormat = new DecimalFormat("##0.00");
    Border border = BorderFactory.createLoweredBevelBorder();
//    GridBagLayout gridBag = new GridBagLayout();
    JCheckBox checkBoxSlider = new JCheckBox();
    JToggleButton selectedButton = new JToggleButton();
    ButtonGroup buttonGroup;
    StsFloatFieldBean sliderValueBean = new StsFloatFieldBean(this, "value", true, null, false);
//    JLabel labelSliderInc = new JLabel();
    StsFloatFieldBean sliderIncBean = new StsFloatFieldBean(this, "increment", true, "Step", false);
    private Color textColor;
    private transient Vector itemListeners;
    private transient Vector changeListeners;
    private double maximum;
    private double minimum;
    private double value;
    private boolean inverted = false;
    private double minIncrement;  // smallest possible increment; other increments are multiples of this
    private double increment;  // current value of the increment
//    private int lastValue = -1;
    public boolean isAdjusting = false;
    private boolean isInitializing = false;
    private boolean showCheckbox = false;

    private Icon selectedIcon = null;
    private Icon deselectedIcon = null;

    static final boolean debug = false;

    public StsSliderBean()
    {
        this(true);
    }

    public StsSliderBean(boolean showCheckbox)
    {
        this(false, showCheckbox);
    }

    // JKF Added this constructor to give the option of having insets - 24MAY2006
    public StsSliderBean(boolean addInsets, boolean showCheckbox)
    {
        super(addInsets);
        this.showCheckbox = showCheckbox;

        try
        {
            initialLayout();
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

    private void initialLayout() throws Exception
    {
        Font font = new java.awt.Font("Dialog", 1, 10);

        checkBoxSlider.setToolTipText("slider on/off");
        checkBoxSlider.setFont(font);
        checkBoxSlider.setMargin(new Insets(0, 0, 0, 0));
//		checkBoxSlider.setPreferredSize(new Dimension(15, 15));
        checkBoxSlider.addItemListener(this);
        setValueLabel("XX"); // sets selectedButton label
        selectedButton.setFont(font);
//        selectedButton.setAlignmentX( (float) 0.5);
        selectedButton.setBorderPainted(false);
//		selectedButton.setPreferredSize(new Dimension(20, 10));
        selectedButton.setMargin(new Insets(0, 0, 0, 0));
        selectedButton.addItemListener(this);

        sliderValueBean.setToolTipText("current slider value");
        sliderValueBean.setBeanFont(font);
//		sliderValueBean.format = valueFormat;
        sliderValueBean.setColumns(6);

        sliderIncBean.setToolTipText("slider value jump");
        sliderIncBean.setBeanFont(font);
//		sliderIncBean.format = incFormat;
        sliderIncBean.setColumns(3);

        slider.setPaintTicks(true);
        slider.setPaintLabels(false);
        slider.addChangeListener(this);
        //Dimension sliderSize = new Dimension(100, 15);
        //slider.setMinimumSize(sliderSize);

        setBorder(border);
        layoutAll();
    }

    private void layoutAll()
    {
        gbc.fill = gbc.NONE;
        gbc.anchor = gbc.WEST;
        gbc.gridwidth = 1;
        if(showCheckbox) addToRow(checkBoxSlider);
        addToRow(selectedButton);
        gbc.gridwidth = 4;
        gbc.fill = gbc.HORIZONTAL;
        addToRow(sliderValueBean);
        gbc.gridwidth = 2;
        addEndRow(sliderIncBean);
        if(showCheckbox)
            gbc.gridwidth = 8;
        else
            gbc.gridwidth = 7;
        add(slider);
        gbc.gridwidth = 1;
    }

    public void showCheckbox(boolean show)
    {
        checkBoxSlider.setEnabled(show);
    }

    public void addToGroup(ButtonGroup buttonGroup)
    {
        this.buttonGroup = buttonGroup;
        buttonGroup.add(selectedButton);
        selectedButton.setBorderPainted(true);
        selectedButton.setSelected(true);
    }

// Event Handling

    /**
     * Sets slider to its range and smallest possible increment.
     *
     * @param min range minimum
     * @param max range maximum
     * @param inc Smallest possible increment; other increments are larger
     *   multiples of this.
     * @param v initial slider setting.
     * @author
     */
    public void initSliderValues(float min, float max, float inc, float v)
    {
        minimum = min;
        maximum = max;
        this.value = v;

        if (minimum >= 1.e30f && maximum >= 1.e30f)
        {
            minimum = 0.f;
            maximum = 100.f;
        }
        inverted = minimum > maximum;

        if(inc == 0.0f) inc = (float)(maximum - minimum);
        minIncrement = inc;
        increment = minIncrement;
        if(inverted)
        {
            if(value > minimum) value = minimum;
            else if(value < maximum) value = maximum;
        }
        else
        {
             if(value < minimum) value = (float)minimum;
             else if(value > maximum) value = (float)maximum;
        }

        Runnable runnable = new Runnable()
        {
            public void run()
            {
                slider.setInverted(inverted);
                isInitializing = true;
                sliderValueBean.setValueAndRange(value, minimum, maximum);
                sliderMax = (int)Math.round((maximum - minimum)/minIncrement);
                int intValue = floatToIntValue((float)value);
                slider.initSliderValues(0, sliderMax, 1, intValue);
                setSliderInc((float)increment);
                isInitializing = false;
            }
        };
        StsToolkit.runLaterOnEventThread(runnable);
    }

    private void setSliderInc(float inc)
    {
        if(inverted) inc = -inc;
        sliderIncBean.setValue(inc);
    }
/*
    public void setLabelFormats(String valueLabel, String incLabel)
    {
        valueFormat = new DecimalFormat(valueLabel);
        incFormat = new DecimalFormat(incLabel);
    }
*/
    private int floatToIntValue(float value)
    {
        return (int)Math.round( (value - minimum)/minIncrement );
    }

    private float intToFloatValue(int value)
    {
        return (float)(minimum + minIncrement*value);
    }

// JavaBeans getters/setters

    public float getValue()
    {
        if(slider == null) return 0.0f;
        return intToFloatValue(slider.getValue());
//        double value = intToFloatValue(slider.getValue());
//        return (float)StsMath.minMax(value, minimum, maximum);
    }

    public void setValue(float value)
    {
        // set the slider
        final int newValue = floatToIntValue(value);
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    slider.setValue(newValue);
                }
            }
        );
    }

    public void setSliderValue(float value)
    {
        sliderValueBean.setValue(value);
    }
    public void incrementSliderValue()
    {
		slider.increment();
        sliderValueBean.increment();
    }
    public void decrementSliderValue()
    {
		slider.decrement();
        sliderValueBean.decrement();
    }

    public void setMaximum(float newMaximum)
    {
        maximum = newMaximum;
    }

    public float getMaximum()
    {
        return (float)maximum;
    }

    public void setMinimum(float newMinimum)
    {
        minimum = newMinimum;
    }

    public float getMinimum()
    {
        return (float)minimum;
    }
/*
    public void setMajorTickSpacing(int newMajorTickSpacing)
    {
        slider.setMajorTickSpacing(newMajorTickSpacing);
    }

    public void setMinorTickSpacing(int newMinorTickSpacing)
    {
        slider.setIncrement(newMinorTickSpacing);
        textSliderInc.setText(String.valueOf(labelFormat.format(
            newMinorTickSpacing)));
    }
*/
    public void setIncrement(float incF)
    {
        if(inverted) incF = -incF;
        if((float)increment == incF) return;

        int inc = (int)Math.round(incF/minIncrement);
        if(inc == 0)
        {
            increment = minIncrement;
            inc = 1;
            sliderIncBean.setValue(minIncrement);
        }
        else
            increment = inc*minIncrement;

        boolean drawTicks = (incF != 0.0f && (maximum - minimum)/incF < 200);
        slider.setIncrement(inc, drawTicks);
        slider.computeNewValue();
    }

    public void setIncrementValue(float incF)
    {
        sliderIncBean.setValue(incF);
    }

    public float getIncrement()
    {
        return (float)increment;
    }

    public void setModel(BoundedRangeModel newModel)
    {
        slider.setModel(newModel);
    }

    public BoundedRangeModel getModel()
    {
        return slider.getModel();
    }

    public void setTextColor(Color newTextColor)
    {
        textColor = newTextColor;
        sliderIncBean.getLabel().setForeground(newTextColor);
        checkBoxSlider.setForeground(newTextColor);
        selectedButton.setForeground(newTextColor);
        selectedButton.setOpaque(true);
    }

    public Color getTextColor()
    {
        return textColor;
    }

    public void setPaintLabels(boolean newPaintLabels)
    {
        slider.setPaintLabels(newPaintLabels);
    }

    public void setPaintTicks(boolean newPaintTicks)
    {
        slider.setPaintTicks(newPaintTicks);
    }

    public void setSelected(boolean enable)
    {
        checkBoxSlider.setSelected(enable);
    }

    public void setCheckBoxModelSelected(boolean enable)
    {
        checkBoxSlider.getModel().setSelected(enable);
    }

    public boolean isSelected()
    {
        return checkBoxSlider.isSelected();
    }

    public boolean isDraggingSlider()
    {
        return slider.isDraggingSlider();
    }

    public boolean isInitializing()
    {
        return isInitializing;
    }

    public void sliderOnly(boolean state)
    {
        reinitializeLayout();

        if (state)
        {
            removeAll();
            gbc.fill = gbc.HORIZONTAL;
            add(slider);
        }
        else
        {
            removeAll();
            layoutAll();
        }
        validate();
        return;
    }

    public JCheckBox getCheckBoxSlider()
    {
        return checkBoxSlider;
    }

    public JLabel getLabelSliderInc()
    {
        return sliderIncBean.getLabel();
    }

    public String getIncrementLabel()
    {
        return sliderIncBean.getLabel().getText();
    }

    public void setIncrementLabel(String label)
    {
        sliderIncBean.getLabel().setText(label);
    }

    public String getValueLabel()
    {
        return selectedButton.getText();
    }

    public void setValueLabel(String label)
    {
        selectedButton.setText(label);
        if(debug) slider.setDebug(label);
    }

    public void setIcons(String selectedName, String deselectedName)
    {
        selectedIcon = StsIcon.createIcon(selectedName + ".gif");
        deselectedIcon = StsIcon.createIcon(deselectedName + ".gif");
    }

// Event listeners

    public synchronized void removeItemListener(ItemListener itemListener)
    {
        if (itemListeners != null && itemListeners.contains(itemListener))
        {
            Vector v = (Vector) itemListeners.clone();
            v.removeElement(itemListener);
            itemListeners = v;
        }
    }

    public synchronized void addItemListener(ItemListener itemListener)
    {
        Vector v = itemListeners == null ? new Vector(2) :
            (Vector) itemListeners.clone();
        if (!v.contains(itemListener))
        {
            v.addElement(itemListener);
            itemListeners = v;
        }
    }

    protected void fireItemStateChanged(ItemEvent e)
    {
        if (itemListeners != null)
        {
            Vector listeners = itemListeners;
            int count = listeners.size();
            for (int i = 0; i < count; i++)
            {
                ( (ItemListener) listeners.elementAt(i)).itemStateChanged(e);
            }
        }
    }

    public synchronized void removeChangeListener(ChangeListener l)
    {
        if (changeListeners != null && changeListeners.contains(l))
        {
            Vector v = (Vector) changeListeners.clone();
            v.removeElement(l);
            changeListeners = v;
        }
    }

    public synchronized void addChangeListener(ChangeListener l)
    {
        Vector v = changeListeners == null ? new Vector(2) :
            (Vector) changeListeners.clone();
        if (!v.contains(l))
        {
            v.addElement(l);
            changeListeners = v;
        }
    }

    protected void fireStateChanged(ChangeEvent e)
    {
        if (changeListeners != null)
        {
            Vector listeners = changeListeners;
            int count = listeners.size();
            for (int i = 0; i < count; i++)
            {
                ( (ChangeListener) listeners.elementAt(i)).stateChanged(e);
            }
        }
    }

    public boolean isPaintLabels()
    {
        return slider.getPaintLabels();
    }

    public boolean isPaintTicks()
    {
        return slider.getPaintTicks();
    }

    public int getMinorTickSpacing()
    {
        return slider.getMinorTickSpacing();
    }

    public int getMajorTickSpacing()
    {
        return slider.getMajorTickSpacing();
    }

    public void itemStateChanged(ItemEvent e)
    {
        e.setSource(this);
        fireItemStateChanged(e);
    }

    public boolean isAdjusting()
    {
		return isAdjusting;
    }

    public boolean adjustingChanged()
    {
        return isAdjusting != slider.getValueIsAdjusting();
    }

    public void stateChanged(ChangeEvent e)
    {
        slider.computeNewValue();

        if (!isInitializing && !checkBoxSlider.isSelected())
        {
            setSelected(true);
        }
		sliderValueBean.setValue(getValue());
//		lastValue = slider.getValue();
		isAdjusting = slider.getValueIsAdjusting();
        if(debug) StsException.systemDebug(this, "stateChanged", "SliderBean " + selectedButton.getText() + " state changed.\n Value: " + slider.getValue() + " Adjusting: " + isAdjusting);
		fireStateChanged(new ChangeEvent(this));
    }

    public void repaint()
    {
        super.repaint();
    }

    public void repaint(Rectangle r)
    {
        validate();
        super.repaint(r);
    }

    public void lock() { slider.lock(); }
    public void unlock() { slider.unlock(); }

    public void highlightSlider(boolean highlight)
    {
        this.putClientProperty("JSlider.isFilled", new Boolean(highlight));
    }

    public void setEditable(boolean enabled)
    {
        //checkBoxSlider.setEnabled(enabled);
        //textSliderValue.setEnabled(enabled);
//        labelSliderInc.setEnabled(enabled);
//        labelSliderInc.setVisible(enabled);
        sliderIncBean.setEditable(enabled);
        sliderIncBean.setVisible(enabled);
        slider.setEnabled(enabled);
        slider.setVisible(enabled);
    }

    public static void main(String[] args)
    {
        try
        {
            UIManager.setLookAndFeel(new com.sun.java.swing.plaf.windows.
                                     WindowsLookAndFeel());
            //    UIManager.setLookAndFeel(new com.sun.java.swing.plaf.motif.MotifLookAndFeel());
            // UIManager.setLookAndFeel(new javax.swing.plaf.metal.MetalLookAndFeel());
        }
        catch (Exception e)
        {
        }

        JFrame frame = new JFrame();

        frame.setSize(new Dimension(480, 100));
        frame.setTitle("S2S 3D Cursor Tool");

        StsSliderBean sliderBean = new StsSliderBean();
//        slider.setLabelFormats("0.000", "0.000");
        sliderBean.initSliderValues(-20.0f, 20.0f, 1.0f, 0.0f);
//        slider.initSliderValues(0.95f, 1.0f, 0.01f, 0.96f);
        frame.getContentPane().add(sliderBean, BorderLayout.CENTER);
        frame.pack();
        frame.setVisible(true);
    }

    public JToggleButton getSelectedButton()
    {
        return selectedButton;
    }

    public void setButtonSelected(boolean b)
    {
        selectedButton.setSelected(b);
        if(b)
        {
            if(selectedIcon != null)
            {
                selectedButton.setIcon(selectedIcon);
            }
            else
            {
                selectedButton.setBackground(textColor);
            //    selectedButton.setForeground(Color.LIGHT_GRAY);
            }
        }
        else
        {
            if(deselectedIcon != null)
            {
                selectedButton.setIcon(deselectedIcon);
            }
            else
            {
                selectedButton.setBackground(Color.LIGHT_GRAY);
            //    selectedButton.setForeground(textColor);
            }
        }
    }

    public void setCheckBoxSelected(boolean b)
        {
            checkBoxSlider.setSelected(b);
    }
    public void setButtonModelSelected(boolean b)
    {
        selectedButton.getModel().setSelected(b);
    }

}
