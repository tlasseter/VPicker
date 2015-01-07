//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import java.awt.*;
import java.awt.event.*;

public class StsMouse
 {
     static StsMouse instance = null;

    // jbw allow remapping
     public static /*final*/ int LEFT = 1;
     public static /*final*/ int MIDDLE = 2;
     public static /*final*/ int RIGHT = 3;
     public static /*final*/ int NONE = 0;

	 public static int POPUP = MIDDLE;
	 public static int VIEW = RIGHT;

     public static String[] buttonLabels = new String[]{"None", "Left", "Middle", "Right"};

     public static final int CLEARED = 0;
     public static final int PRESSED = 1;
     public static final int DRAGGED = 2;
     public static final int RELEASED = 3;
     public static final int WHEEL = 4;
	 // jbw for more subtle handling of mouse events
	 // public static final int CLICKED = 5;
	 // public static final int DOUBLECLICKED = 6;
     // tjl: event sequence is unreliable, so we need to handle CLICK and DOUBLECLICK with a timer
     public static boolean wasClicked = false;

     static private int[] buttonStates = new int[]{CLEARED, CLEARED, CLEARED, CLEARED, CLEARED, CLEARED}; // initialize all to CLEARED

     public static String[] stateLabels = new String[]{"Cleared", "Pressed", "Dragged", "Released", "Wheel", "Click", "Doubleclcick"};

     public static int currentButton = NONE;

     public static boolean shiftDown = false;
     public static boolean controlDown = false;
     public static boolean altDown = false;
	 public static boolean remapped = false;

     public static StsMousePoint mousePoint = new StsMousePoint();
     public static StsMousePoint mouseDelta = new StsMousePoint();

     public static long pressedTime = 0L;
     public static long clickInterval = 250;

     private static final boolean debug = false;

     private StsMouse()
     {
     }

     static public StsMouse getInstance()
     {
         if (instance == null) instance = new StsMouse();
         return instance;
     }

     //    static private int currentKeyPressed = KeyEvent.CHAR_UNDEFINED;

     /*	public StsMouse()
     {
     }
     */
	public void remap(boolean b)
	{
		remapped = b;
		if (remapped)
		{
			POPUP = RIGHT;
			VIEW = MIDDLE;
		}
		else
		{
			POPUP = MIDDLE;
			VIEW = RIGHT;
		}
	}
     /**
      * This setState is called by a component supporting enableEvents
      * such as StsColorGridFrame with a corresponding processMouseEvent
      * method.  We get the actionID and call setState(e, actionID).
      * If the component supports only a MouseListener interface such as
      * the GLComponent, individual mouse methods must be supported/called
      * such as mousePressed, etc.  In this case, we set the actionID there
      * and call setState(e, actionID) directly.  We could of course call
      * setState(e) and reextract the actionID...
      */
     public void setState(MouseEvent e)
     {
         switch (e.getID())
         {
             // ignore CLICK event:  handled by RELEASE event with timer
             case MouseEvent.MOUSE_CLICKED:
				 // setState(e, CLICKED);
                 return;

             case MouseEvent.MOUSE_PRESSED:
                 setState(e, PRESSED);
                 break;
             case MouseEvent.MOUSE_DRAGGED:
                 setState(e, DRAGGED);
                 break;
             case MouseEvent.MOUSE_RELEASED:
                 setState(e, RELEASED);
                 break;
             case MouseEvent.MOUSE_WHEEL:
                 setState(e, WHEEL);
                 break;

             default:
                 setState(e, CLEARED);
         }
     }

     public void setState(MouseEvent e, int actionID)
     {
         //  When mouse dragged, button ID is not captured, so we have to
         //  maintain state of current button.

         // Button might be PRESSED, DRAGGED, and RELEASED or
         // it might be PRESSED, RELEASED, and CLICKED.
         // So NEVER call this routine if CLICKED. We will use RELEASED as
         // our trigger event.

         shiftDown = e.isShiftDown();
         controlDown = e.isControlDown();
         altDown = e.isAltDown();

         int mods = e.getModifiers();

         if (mods != 0)
         {
             //            if (shiftDown && ((mods & InputEvent.BUTTON1_MASK) != 0))
             //                currentButton = MIDDLE;   // shift-left btn same as middle btn
             if (shiftDown && ((mods & InputEvent.BUTTON3_MASK)) != 0)
                 currentButton = MIDDLE;   // shift-right btn same as middle btn
             else if ((mods & InputEvent.BUTTON1_MASK) != 0)
                 currentButton = LEFT;
             else if ((mods & InputEvent.BUTTON3_MASK) != 0)
                 currentButton = RIGHT;
             else if ((mods & InputEvent.BUTTON2_MASK) != 0)
                 currentButton = MIDDLE;   // never gets here w/ 2-btn mouse
             else
             {
                 currentButton = NONE;
                 return;
             }
         }

         if(actionID == PRESSED)
         {
            pressedTime = System.currentTimeMillis();
            wasClicked = false;
         }
         else if(actionID == RELEASED)
         {
             long releasedTime = System.currentTimeMillis();
             long interval = releasedTime - pressedTime;
             wasClicked = interval < clickInterval;
         }

         if(debug) System.out.println("Button " + buttonLabels[currentButton] + " set to " + stateLabels[actionID]);
         buttonStates[currentButton] = actionID;
         setMouseCoordinates(e, actionID);
     }

     static public boolean isShiftDown()
     {
         return shiftDown;
     }

     static public boolean isAltDown()
     {
         return altDown;
     }

     static public boolean isControlDown()
     {
         return controlDown;
     }

     static public boolean isLeftButtonDragging()
     {
         return buttonStates[LEFT] == DRAGGED;
     }

     static public boolean isLeftButtonDown()
     {
//         System.out.println("left button state: " + stateLabels[buttonStates[LEFT]]);
         return buttonStates[LEFT] == PRESSED;
     }

     static public boolean isRightButtonDragging()
     {
         return buttonStates[RIGHT] == DRAGGED;
     }

     static public boolean isRightButtonDown()
     {
         // System.out.println("right button state: " + stateLabels[buttonStates[RIGHT]]);
         return buttonStates[RIGHT] == PRESSED;
     }

     static public boolean isButtonStateReleasedOrClicked(int buttonState)
     {
         return buttonState == RELEASED;
     }

     static public boolean isButtonStateClicked(int buttonState)
     {
         return buttonState == RELEASED && wasClicked;
     }
     /*
         static public void keyPressed(int keyChar)
         {
             if(keyChar == KeyEvent.VK_Z)
                 currentKeyPressed = keyChar;
             else if(keyChar == KeyEvent.VK_X)
                 currentKeyPressed = keyChar;
         }

         static public void keyReleased(int keyChar)
         {
             currentKeyPressed = KeyEvent.CHAR_UNDEFINED;
         }

         static public int getKeyPressed()
         {
             return currentKeyPressed;
         }
     */
     public void setMouseCoordinates(MouseEvent e, int actionID)
     {
         if (actionID == DRAGGED)
         {
             int newX = e.getX();
             int newY = e.getY();

             mouseDelta.x = newX - mousePoint.x;
             mouseDelta.y = -newY + mousePoint.y;

             mousePoint.x = newX;
             mousePoint.y = newY;
         }
         else if (actionID == WHEEL)
         {
             mouseDelta.y = -5 * ((MouseWheelEvent) e).getWheelRotation();
             mouseDelta.x = 0;
         }
         else
         {
             mouseDelta.x = 0;
             mouseDelta.y = 0;
             mousePoint.x = e.getX();
             mousePoint.y = e.getY();
         }
     }

     public void printState()
     {
         System.out.println(toString());
     }

     public String toString()
     {
         return new String(buttonLabels[currentButton] + " " + stateLabels[buttonStates[currentButton]]);
     }

     public StsMousePoint getMousePoint()
     {
         return mousePoint;
     }

     public StsMousePoint getMousePointGL(Rectangle window)
     {
         return mousePoint.lowerLeftCoor(window);
     }

     public StsMousePoint getMousePointGL(int[] viewPort)
     {
         return mousePoint.lowerLeftCoor(viewPort);
     }

     public int getX() { return mousePoint.x; }
     public int getY() { return mousePoint.y; }

     public StsMousePoint getMouseDelta()
     {
         return mouseDelta;
     }

     public boolean isButtonStatePressed(int buttonNumber)
     {
         if(buttonStates[buttonNumber] != PRESSED) return false;
         return true;
     }

     public boolean isButtonStateDragging(int buttonNumber)
     {
         if(buttonStates[buttonNumber] != DRAGGED) return false;
         return true;
     }

     public boolean isButtonStatePressedOrDragged(int buttonNumber)
     {
         if(buttonStates[buttonNumber] != PRESSED) return false;
         if(buttonStates[buttonNumber] != DRAGGED) return false;
         return true;
     }

     public boolean isButtonStateReleased(int buttonNumber)
     {
         if(buttonStates[buttonNumber] != RELEASED) return false;
         clearButtonState(buttonNumber);
         return true;
     }

     public int getButtonStateCheckClear(int buttonNumber)
     {
         int buttonState = buttonStates[buttonNumber];
         if(buttonState == RELEASED)
            clearButtonState(buttonNumber);
         return buttonState;
     }

     public int getButtonState(int buttonNumber)
     {
         int buttonState = buttonStates[buttonNumber];
         return buttonState;
     }

     public void clearButtonState(int buttonNumber)
     {
         buttonStates[buttonNumber] = CLEARED;
     }

     public int getCurrentButton()
     {
         return currentButton;
     }

     public boolean isButtonDown()
     {
         if (currentButton == NONE) return false;
         return isButtonDown(currentButton);
     }

     public boolean isButtonDown(int button)
     {
         int buttonState = buttonStates[button];
         return buttonState == PRESSED || buttonState == DRAGGED;
     }
 }
