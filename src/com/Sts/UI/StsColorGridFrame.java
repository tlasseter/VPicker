
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;

public class StsColorGridFrame extends JFrame implements MouseListener, MouseMotionListener
{
	static int buttonWidth = 50;
    static int buttonHeight = 50;

    protected int nRows, nCols;
    protected StsSpectrum spectrum = null;
   	protected int nButtons;
    protected JPanel[] buttons = null;
    protected int currentButtonNo;

    static BevelBorder raised = new BevelBorder(BevelBorder.RAISED);
    static BevelBorder lowered = new BevelBorder(BevelBorder.LOWERED);

    static private StsMouse mouse = StsMouse.getInstance();

    /** Construct a JFrame with a specified number of rows and columns of
      * JPanel buttons with the colors from the spectrum and the current Button
      * selected.
      * @param nRows number of rows of buttons
      * @param nCols number of columns of buttons
      * @spectrum a spectrum of StsColors (nRows*nCols)
      * @currentButton index of current button to be selected
      *
      */
	public StsColorGridFrame(int nRows, int nCols, StsSpectrum spectrum, int currentButtonNo)
	{
       /** We want mouse and window events */
       /*
		enableEvents(AWTEvent.WINDOW_EVENT_MASK |
                     AWTEvent.MOUSE_EVENT_MASK |
                     AWTEvent.MOUSE_MOTION_EVENT_MASK);

        */

        addMouseListener(this);

    	this.nRows = nRows;
        this.nCols = nCols;
        this.spectrum = spectrum;
        this.currentButtonNo = currentButtonNo;

    	nButtons = nRows*nCols;

        StsColor[] colors = spectrum.getStsColors();
        int nColors = colors.length;
    /*
        if(colors.length != nButtons)
        {
            System.err.println("Spectrum size is: " + colors.length +
                               ", but there are " + nButtons + " buttons provided.");
            return;
        }
    */
		try
		{
		    this.setTitle("Select Current Color");
            this.setSize(new Dimension(nCols*buttonWidth, nRows*buttonHeight));

            Container pane = getContentPane();
  	        GridLayout gridLayout = new GridLayout();
            pane.setLayout(gridLayout);
		    gridLayout.setRows(nRows);
		    gridLayout.setColumns(nCols);

  		    buttons = new JPanel[nButtons];

            /** Add buttons with color, raised bevel, and number label */
            for(int n = 0; n < nColors; n++)
            {
        	    buttons[n] = new JPanel();
			    pane.add(buttons[n], null);

                buttons[n].setBackground(colors[n].getColor());
                buttons[n].setBorder(raised);

                if(nColors <= 32)
                {
                    JLabel label = new JLabel(Integer.toString(n));
                    label.setForeground(Color.black);
                    buttons[n].add(label, "Center");
                }
                // Would be nice to use toolTips, but they swallow
                // mouse events.  Need to work on it.  One issue is
                // that we need to go to full swing compliance if we
                // want ToolTips in the Win3d window, but Swing/JDK1.2
                // doesn't work with Magician: lightweight/heavyweight
                // problems.
                buttons[n].setToolTipText(colors[n].getHSBString());
            }

            /** Reset selected button with lowered bevel */
            buttons[currentButtonNo].setBorder(lowered);

            StsToolkit.centerComponentOnScreen(this);
		    setVisible(true);
 		}
		catch (Exception e)
		{
			StsException.outputException("Failed building StsColorGridFrame.", e, StsException.WARNING);
		}
	}

// METHODs: MouseListener interface: Clicked, Entered, Exited, Pressed, Released
// Magician currently doesn't support enableEvent design: need to you Listener for now...

    public void mouseClicked(MouseEvent e)
	{
    // Ignore Clicked since we always get Released

    //    performMouseAction(e, MouseEvent.MOUSE_CLICKED);
    }

    public void mousePressed(MouseEvent e)
	{
    	try
        {
        	performMouseAction(e, StsMouse.PRESSED);
        }
        catch (StsException ex)
        {
        	ex.warning();
            return;
        }
    }

    public void mouseReleased(MouseEvent e)
	{
     	try
        {
        	performMouseAction(e, StsMouse.RELEASED);
        }
         catch (StsException ex)
        {
        	ex.warning();
            return;
        }
    }


    public void mouseDragged(MouseEvent e)
	{
     	try
        {
        	performMouseAction(e, StsMouse.DRAGGED);
        }
         catch (StsException ex)
        {
        	ex.warning();
            return;
        }
    }

    private void performMouseAction(MouseEvent e, int mouseActionID) throws StsException
    {
        mouse.setState(e, mouseActionID);

        if(mouse.getCurrentButton() == StsMouse.LEFT)
        {
            int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

            if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
                setButton(mouse);
        }

        if(mouse.getCurrentButton() == StsMouse.VIEW)
        {
            if(mouse.getButtonStateCheckClear(StsMouse.VIEW) == StsMouse.PRESSED)
                selectNewColor(mouse);
        }
    }

	public void mouseMoved(MouseEvent e)
	{
	}

    public void mouseEntered(MouseEvent e)
	{
	}

    public void mouseExited(MouseEvent e)
	{
	}

/*
    public void processMouseEvent(MouseEvent e)
    {
        handleMouseEvent(e);
    }

    public void processMouseMotionEvent(MouseEvent e)
    {
       handleMouseEvent(e);
    }

    private void handleMouseEvent(MouseEvent e)
    {
        mouse.setState(e);

        if(mouse.getCurrentButton() == StsMouse.LEFT)
        {
            int leftButtonState = mouse.getButtonState(StsMouse.LEFT);

            if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
                setButton(mouse);
            else
                mouse.clearButtonState(StsMouse.LEFT);
        }

        if(mouse.getCurrentButton() == StsMouse.VIEW)
        {
            if(mouse.getButtonState(StsMouse.VIEW) == StsMouse.PRESSED)
                selectNewColor(mouse);
        }
    }
*/
    /** This is very hokey: a mouseEvent in the top bar causes a mouseEvent
      * here.  We need to compute the geometry of the buttons and bar to
      * filter this out.  Would be nicer if mouseEvents were restricted to
      * the inside region of the window.  It should be easier to get the
      * border size than I have done here.
      */
    private void setButton(StsMouse mouse)
    {
        int buttonNo = getButtonNo(mouse);
        if(buttonNo < 0) return;

        buttonSelected(buttonNo);
        spectrum.setCurrentColorIndex(buttonNo);
    }

    private int getButtonNo(StsMouse mouse)
    {
        Point windowPoint = this.getLocation();
        Point buttonPoint = buttons[0].getLocationOnScreen();
        int leftBorder = buttonPoint.x - windowPoint.x;
        int topBorder = buttonPoint.y - windowPoint.y;

        StsMousePoint mousePoint = mouse.getMousePoint();
        Rectangle buttonRectangle = buttons[0].getBounds();

        if(mousePoint.y < topBorder) return -1;

        int buttonCol = (mousePoint.x - leftBorder)/buttonRectangle.width;
        int buttonRow = (mousePoint.y - topBorder)/buttonRectangle.height;

        return buttonRow*nCols + buttonCol;
    }

    /** Button has been selected for this spectrum:
      * set current color index
      */
    public void buttonSelected(int index)
    {
        buttons[currentButtonNo].setBorder(raised);
        buttons[currentButtonNo].repaint();

        buttons[index].setBorder(lowered);
        buttons[index].repaint();

        currentButtonNo = index;
    }

    private void selectNewColor(StsMouse mouse)
    {
        int buttonNo = getButtonNo(mouse);
        if(buttonNo < 0) return;

        StsColor[] colors = spectrum.getStsColors();
        StsColor buttonColor = colors[buttonNo];

        Color currentColor = new Color(buttonColor.red, buttonColor.green, buttonColor.blue);

        Color newColor = JColorChooser.showDialog(buttons[buttonNo], "Select a new color", currentColor);
        if(newColor == null) return;

        colors[buttonNo] = new StsColor(newColor);
        buttons[buttonNo].setBackground(newColor);
    }

//Overriden so we can exit on System Close

	protected void processWindowEvent(WindowEvent e)
	{
		super.processWindowEvent(e);
		if (e.getID() == WindowEvent.WINDOW_CLOSING)
        	setVisible(false);
	}
}


