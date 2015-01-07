
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

public class StsCursor3dDialog extends JFrame
	implements ActionListener, ItemListener, ChangeListener
{
	public static final String closeAction = "close";
	public static final String selectSliderAction = "selectSlider";

    //Construct the frame
    private JPanel panelSliders = new JPanel();
    private JPanel panelButtons = new JPanel();

    private StsSliderBean sliderX = new StsSliderBean();
    private StsSliderBean sliderY = new StsSliderBean();
    private StsSliderBean sliderZ = new StsSliderBean();

    private JPanel jPanel2 = new JPanel();
    private JPanel jPanel1 = new JPanel();
    private JPanel jPanel3 = new JPanel();

    private JRadioButton radioButtonX = new JRadioButton();
    private JRadioButton radioButtonY = new JRadioButton();
    private JRadioButton radioButtonZ = new JRadioButton();

    private GridLayout gridLayout1 = new GridLayout();
    private GridLayout gridLayout2 = new GridLayout();
    private JRootPane jRootPane1 = new JRootPane();
    private BorderLayout borderLayout1 = new BorderLayout();
    private FlowLayout flowLayout1 = new FlowLayout();
    private JPanel jPanel4 = new JPanel();
    private BorderLayout verticalFlowLayout1 = new BorderLayout();
    private BorderLayout verticalFlowLayout2 = new BorderLayout();

    private JLabel labelCurrentPlane = new JLabel();
    private JButton buttonClose = new JButton();
    private JButton buttonHelp = new JButton();

    private StsSliderBean selectedSlider;
    private GridLayout gridLayout3 = new GridLayout();
    private BorderLayout borderLayout2 = new BorderLayout();

	private transient Vector actionListeners = null;
	private transient Vector itemListeners = null;
	private transient Vector changeListeners = null;

    private boolean isInitializing = false;

   	public StsCursor3dDialog()
    {
        try
        {
            jbInit();
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

//Component initialization

    private void jbInit() throws Exception
    {
        this.getContentPane().setLayout(borderLayout2);
        this.setSize(new Dimension(418, 263));
        panelButtons.setLayout(gridLayout2);

        sliderX.setTextColor(Color.red);
        sliderY.setTextColor(Color.green);
        sliderZ.setTextColor(Color.blue);

        sliderX.getCheckBoxSlider().setText("X-Slider Value:");
        sliderY.getCheckBoxSlider().setText("Y-Slider Value:");
        sliderZ.getCheckBoxSlider().setText("Z-Slider Value:");

        sliderX.getCheckBoxSlider().setSelected(false);
        sliderY.getCheckBoxSlider().setSelected(false);
        sliderZ.getCheckBoxSlider().setSelected(false);

        sliderX.getLabelSliderInc().setText("X-Increment:");
        sliderY.getLabelSliderInc().setText("Y-Increment:");
        sliderZ.getLabelSliderInc().setText("Z-Increment:");

        jPanel2.setLayout(flowLayout1);
        radioButtonX.setMaximumSize(new Dimension(53, 23));
        radioButtonX.setText("X");
        radioButtonX.setFont(new Font("Dialog", 1, 12));

        radioButtonY.setMaximumSize(new Dimension(53, 23));
        radioButtonY.setText("Y");
        radioButtonY.setPreferredSize(new Dimension(33, 23));
        radioButtonY.setFont(new Font("Dialog", 1, 12));
        radioButtonY.setMinimumSize(new Dimension(33, 23));

        radioButtonZ.setMaximumSize(new Dimension(53, 23));
        radioButtonZ.setText("Z");
        radioButtonZ.setPreferredSize(new Dimension(33, 23));
        radioButtonZ.setFont(new Font("Dialog", 1, 12));
        radioButtonZ.setMinimumSize(new Dimension(33, 23));

        labelCurrentPlane.setPreferredSize(new Dimension(131, 15));
        labelCurrentPlane.setText("Current Plane Selected");
        labelCurrentPlane.setFont(new Font("Dialog", 1, 12));
        labelCurrentPlane.setAlignmentX((float) 0.5);
        jPanel4.setLayout(verticalFlowLayout1);
//        verticalFlowLayout1.setAlignment(1);
//        verticalFlowLayout2.setAlignment(1);
        buttonClose.setText("Close");
        buttonClose.setPreferredSize(new Dimension(100, 23));
        buttonClose.setFont(new Font("Dialog", 1, 12));

        buttonHelp.setText("Help");
        buttonHelp.setFont(new Font("Dialog", 1, 12));
        gridLayout3.setRows(3);
        gridLayout3.setColumns(1);
        jRootPane1.getContentPane().setLayout(verticalFlowLayout2);
        jPanel3.setLayout(gridLayout1);
        jPanel1.setLayout(borderLayout1);

        panelSliders.setLayout(gridLayout3);
        this.setTitle("S2S 3D Cursor Tool");
        this.getContentPane().add(panelSliders, BorderLayout.CENTER);

        panelSliders.add(sliderX, null);
        panelSliders.add(sliderY, null);
        panelSliders.add(sliderZ, null);

        this.getContentPane().add(panelButtons, BorderLayout.SOUTH);
        panelButtons.add(jPanel4, null);
        jPanel4.add(buttonClose, BorderLayout.CENTER);
        panelButtons.add(jPanel2, null);
        jPanel2.add(jPanel1, null);
        jPanel1.add(jPanel3, BorderLayout.NORTH);
        jPanel3.add(radioButtonX, null);
        jPanel3.add(radioButtonY, null);
        jPanel3.add(radioButtonZ, null);
        jPanel1.add(labelCurrentPlane, BorderLayout.SOUTH);
        panelButtons.add(jRootPane1, null);
        jRootPane1.getContentPane().add(buttonHelp, BorderLayout.CENTER);

        ButtonGroup buttonGroupXYZ = new ButtonGroup();
		gridLayout3.setVgap(2);
		gridLayout3.setHgap(2);
        buttonGroupXYZ.add(radioButtonX);
        buttonGroupXYZ.add(radioButtonY);
        buttonGroupXYZ.add(radioButtonZ);


// Add listeners

        buttonClose.addActionListener(this);
        buttonClose.setActionCommand(closeAction);

        radioButtonX.addActionListener(this);
        radioButtonX.setActionCommand(selectSliderAction);

        radioButtonY.addActionListener(this);
        radioButtonY.setActionCommand(selectSliderAction);

        radioButtonZ.addActionListener(this);
        radioButtonZ.setActionCommand(selectSliderAction);

        sliderX.addChangeListener(this);
        sliderX.addItemListener(this);

        sliderY.addChangeListener(this);
        sliderY.addItemListener(this);

        sliderZ.addChangeListener(this);
        sliderZ.addItemListener(this);

        pack();
   }

    public void initValues(StsCursor3d cursor3d)
    {
    	isInitializing = true;
        StsRotatedGridBoundingBox rotatedBoundingBox = cursor3d.rotatedBoundingBox;

        initSliderValues(sliderX, (int)rotatedBoundingBox.xMin, (int)rotatedBoundingBox.xMax, (int)rotatedBoundingBox.xInc);
        initSliderValues(sliderY, (int)rotatedBoundingBox.yMin, (int)rotatedBoundingBox.yMax, (int)rotatedBoundingBox.yInc);
        initSliderValues(sliderZ, (int)rotatedBoundingBox.zMin, (int)rotatedBoundingBox.zMax, (int)rotatedBoundingBox.zInc);
    	isInitializing = false;
    }

    private void initSliderValues(StsSliderBean slider, int min, int max, int increment)
    {
        int value = (max + min)/2;
        slider.initSliderValues(min, max, increment, value);
    }

    public StsSliderBean getSliderX() { return sliderX; }
    public StsSliderBean getSliderY() { return sliderY; }
    public StsSliderBean getSliderZ() { return sliderZ; }

    void selectSlider(StsSliderBean slider)
    {
        slider.getCheckBoxSlider().setSelected(true);
    }

	public synchronized void removeActionListener(ActionListener l)
	{
		if (actionListeners != null && actionListeners.contains(l))
		{
			Vector v = (Vector) actionListeners.clone();
			v.removeElement(l);
			actionListeners = v;
		}
	}

	public synchronized void addActionListener(ActionListener l)
	{
		Vector v = actionListeners == null ? new Vector(2) : (Vector) actionListeners.clone();
		if (!v.contains(l))
		{
			v.addElement(l);
			actionListeners = v;
		}
	}

	protected void fireActionPerformed(ActionEvent e)
	{
		if (actionListeners != null)
		{
			Vector listeners = actionListeners;
			int count = listeners.size();
			for (int i = 0; i < count; i++)
			  ((ActionListener) listeners.elementAt(i)).actionPerformed(e);
		}
	}

	public synchronized void removeItemListener(ItemListener l)
	{
		if (itemListeners != null && itemListeners.contains(l))
		{
			Vector v = (Vector) itemListeners.clone();
			v.removeElement(l);
			itemListeners = v;
		}
	}

	public synchronized void addItemListener(ItemListener l)
	{
		Vector v = itemListeners == null ? new Vector(2) : (Vector) itemListeners.clone();
		if (!v.contains(l))
		{
			v.addElement(l);
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
			  ((ItemListener) listeners.elementAt(i)).itemStateChanged(e);
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
		Vector v = changeListeners == null ? new Vector(2) : (Vector) changeListeners.clone();
		if (!v.contains(l))
		{
			v.addElement(l);
			changeListeners = v;
		}
        fireStateChanged(new ChangeEvent(this));
	}

	protected void fireStateChanged(ChangeEvent e)
	{
		if (changeListeners != null)
		{
			Vector listeners = changeListeners;
			int count = listeners.size();
			for (int i = 0; i < count; i++)
			  ((ChangeListener) listeners.elementAt(i)).stateChanged(e);
		}
	}

	public void actionPerformed(ActionEvent e)
	{
    	if( e.getActionCommand().equals(closeAction) )
	        setVisible(false);
        else if( e.getActionCommand().equals(selectSliderAction) )
        {
            if( e.getSource() == radioButtonX )
            	selectedSlider = sliderX;
            else if( e.getSource() == radioButtonY )
            	selectedSlider = sliderY;
            else if( e.getSource() == radioButtonZ )
            	selectedSlider = sliderZ;

            if( selectedSlider != null )
	        	selectSlider(selectedSlider);
        }
		fireActionPerformed(e);
	}

	public void itemStateChanged(ItemEvent e)
	{
        if( !(e.getItem() instanceof StsSliderBean) ) return;
        StsSliderBean source = (StsSliderBean) e.getItem();
        int[] dirs = new int[] { StsCursor3d.XDIR, StsCursor3d.YDIR, StsCursor3d.ZDIR };
        StsSliderBean[] sliders = new StsSliderBean[] { sliderX, sliderY, sliderZ };
        JRadioButton[] buttons = new JRadioButton[] { radioButtonX, radioButtonY, radioButtonZ };

		// which direction slider
        int direction = StsParameters.NO_MATCH;
		for( int i=0; i<sliders.length; i++ )
        {
	        if( source == sliders[i] )
            {
                direction = dirs[i];
                break;
            }
        }

        // set the current radio button
       	if( e.getStateChange() == ItemEvent.SELECTED )
		{
            for( int i=0; i<sliders.length; i++ )
            {
                if( source == sliders[i] )
                {
	                buttons[i].setSelected(true);
                    break;
                }
            }
            selectedSlider = source;
        }

        // find a new current
        else if( e.getStateChange() == ItemEvent.DESELECTED )
        {
            for( int i=0; i<dirs.length; i++ )
            {
        		if( dirs[i] == direction ) continue;
                if( sliders[i].getCheckBoxSlider().isSelected() )
                {
					buttons[i].setSelected(true);
                    break;
                }
            }
        }

        ItemEvent event = new ItemEvent(e.getItemSelectable(), direction,
                                        source, e.getStateChange());
        fireItemStateChanged(event);
	}

	public void stateChanged(ChangeEvent e)
	{
        if( !(e.getSource() instanceof StsSliderBean) ) return;
        StsSliderBean source = (StsSliderBean) e.getSource();
        if( !isInitializing )
        {
            selectSlider(source);
            selectedSlider = source;
        }
        fireStateChanged(e);
	}
}


