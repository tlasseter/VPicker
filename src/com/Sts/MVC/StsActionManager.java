package com.Sts.MVC;

import com.Sts.Actions.*;
import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.Workflow.StsWorkflowTreeNode;

import java.lang.reflect.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Dec 14, 2009
 * Time: 3:37:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsActionManager
{
    /** convenience copy of model */
    protected StsModel model;
    /** actionClass currently being executed */
    private Class actionClass;//    private char currentKeyPressed = ' ';
    /**
     * list of actionClasses on stack queue
     * Generally there is only one action, but it may be interrupted
     * while a passive action (surface readout for example) is run.
     */
    private ArrayList actionList = new ArrayList();
    /** name of workflow step being executed */
    private HashMap actionWorkflowHashtable = new HashMap(25);
    /** current workflow step name (assigned to any action created during this workflow step) */
    private String currentWorkflowStepName = "Default";
    /** Do we want actions to wait or create another thread. */
    public boolean waitOnThread = false;
    /** Should thread still be waiting */
    public boolean waiting = false;
    public boolean hasFocus = false;
    static public final boolean debug = false;
    static public final boolean actionDebug = false;
    /** Statics for logging */
    static final int END = 0;
    static final int START = 1;
    static final int ABORT = 2;

    public StsActionManager(StsModel model)
    {
        this.model = model;
    }

    public StsModel getModel() { return model; }

    public void setModel(StsModel model) { this.model = model; }

    public StsAction getCurrentAction()
    {
        int nActions = actionList.size();
        if (nActions == 0) return null;
        return (StsAction) actionList.get(nActions - 1);
    }

    public void viewChanged(Object viewObject)
    {
        int nActions = actionList.size();
        for (int n = 0; n < nActions; n++)
        {
            StsAction action = (StsAction) actionList.get(n);
            action.viewChanged(viewObject);
        }
    }

    /** return true if action doesn't exist or is ended ok */
    public boolean endAllAction()
    {
        StsAction currentAction = getCurrentAction();
        if (currentAction == null) return true;

        // if this is a passive action (canInterrupt = true) terminate it
        if (currentAction.canInterrupt()) return true;

        // if this action is not running an open transaction, don't bother the user
//        if(!currentAction.hasOpenTransaction()) return true;

        //TODO I think this is causing DB problems.  Commented out for now.  Check.  TJL 3/13/09
        /*
        boolean ok = StsMessage.questionValue(model.win3d, currentAction.getClass().getName() +
            " has not been terminated, Do you wish to save transactions?\n" +
            "If not, changes will be lost.");
        if (!ok)
            abortCurrentAction();
        else
            terminateCurrentAction();
        */
        // terminate is probably writing to the DB when it can't be done properly. Fix so this works ok.
        abortCurrentAction();
        return false;
    }

    /** unconditionally terminate currentAction */
    private void terminateCurrentAction()
    {
        StsAction currentAction = getCurrentAction();
        if (currentAction == null) return;
        if (actionDebug) StsException.systemDebug(this, "terminateCurrentAction", "currentAction: " + currentAction.getName());
        if (!currentAction.end()) abortCurrentTransaction();
//        refreshObjectTree();
//		fireChangeEvent();
        logUsage(END, currentAction);

        model.commit();

//        checkEndAuxiliaryWindowAction(currentAction);

        if (currentAction.isRepeatable())
            startAction(actionClass, currentAction);
        else
            removeCurrentAction();
    }

    /** terminate current action and start new one if repeatable */
    public void endCurrentAction()
    {
        StsAction currentAction = getCurrentAction();
        if (currentAction == null) return;
        terminateCurrentAction();
    }

    /** terminate current action and start new one if repeatable */
    public void cancelCurrentAction()
    {
        StsAction currentAction = getCurrentAction();
        if (currentAction == null) return;

        if (currentAction.isRepeatable())
        {
            StsAction saveCurrentAction = currentAction;
            abortCurrentAction();
            startAction(actionClass, saveCurrentAction);
        }
        else
            abortCurrentAction();
    }

    public void addActionWorkflow(Class actionClass)
    {
        getWorkflowStepName(actionClass);
        return;
    }

    private String getWorkflowStepName(StsAction action)
    {
        String workflowStepName = (String) actionWorkflowHashtable.get(action.getClass().getName());
        if (workflowStepName != null) return workflowStepName;
        actionWorkflowHashtable.put(action.getClass().getName(), currentWorkflowStepName);
        return currentWorkflowStepName;
    }

    private String getWorkflowStepName(Class actionClass)
    {
        String actionClassName = actionClass.getName();
        String workflowStepName = (String) actionWorkflowHashtable.get(actionClassName);
        if (workflowStepName != null) return workflowStepName;
        actionWorkflowHashtable.put(actionClassName, currentWorkflowStepName);
        return currentWorkflowStepName;
    }

    private Class[] getArgTypes(Object[] args)
    {
        if (args == null) return new Class[0];
        int nArgs = args.length;
        if (nArgs == 0) return new Class[0];
        Class[] argTypes = new Class[nArgs];
        for (int n = 0; n < nArgs; n++)
            argTypes[n] = args[n].getClass();
        return argTypes;
    }

    private Object[] getArgs(Object[] args)
    {
        if (args == null || args.length == 0)
            return new Object[]{this};
        else
        {
            int nArgs = args.length;
            Object[] newArgs = new Object[1 + nArgs];
            newArgs[0] = this;
            System.arraycopy(args, 0, newArgs, 1, nArgs);
            return newArgs;
        }
    }

    public synchronized void threadWait() throws InterruptedException
    {
//        if(!waitOnThread) return;
        waiting = true;
        String threadName = Thread.currentThread().getName();
//        System.out.println("Thread " + threadName + " is waiting.");
        while (waiting) { wait(); }
    }

    public synchronized void threadComplete() throws InterruptedException
    {
        String threadName = Thread.currentThread().getName();
//        System.out.println("Thread " + threadName + " is complete; notifying waiting thread.");
        waiting = false;
        notifyAll();
    }

    /**
     * Create an action starrted from a workflow step.
     * If the action doesn't already have a workflowStepName, assign it the passed in name.
     * Set the current workflow step name to this name.  Subsequent actions initiated as
     * part of this workflow step will be assigned this workflow step name if they don't
     * already have a workflow step name.
     */
    public StsAction newWorkflowAction(Class actionClass, String workflowStepName)
    {
        try
        {
            if (actionClass == null) return null;
            this.actionClass = actionClass;

            accumulateActions(actionClass, workflowStepName);

            return newAction(actionClass);
        }
        catch (Exception e)
        {
            StsException.outputException("StsActionManager.newWorkflowAction() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    /**
     * Launch a wizard from non-standard (ie. to workflow panel) location
     * @param wizardClass
     * @param wizardName
     * @return
     */
    public boolean launchWizard(String wizardClass, String wizardName)
    {
         ClassLoader classLoader =  StsWorkflowTreeNode.class.getClassLoader();
         try
         {
             Class actionClass = classLoader.loadClass(wizardClass);
             StsWizard wizard = (StsWizard)newWorkflowAction(actionClass, wizardName);
             checkStartAction(wizard);
             return true;
         }
         catch(Exception e)
         {
             StsException.outputException("Unable to launch wizard: " + wizardName, e, StsException.WARNING);
             return false;
         }
    }

    /** create instance of actionClass but don't start it */
    public void accumulateActions(Class actionClass, String workflowStepName)
    {
        String actionClassName = actionClass.getName();
        String actionWorkflowStepName = (String) actionWorkflowHashtable.get(actionClassName);
        if (actionWorkflowStepName == null)
        {
            actionWorkflowHashtable.put(actionClassName, workflowStepName);
            currentWorkflowStepName = workflowStepName;
        }

        return;
    }

    /** create instance of actionClass but don't start it */
    public StsAction newAction(Class actionClass)
    {
        if (actionClass == null) return null;
        this.actionClass = actionClass;

        accumulateActions(actionClass, "None");

        return newAction(actionClass, null);
    }

    /** create instance of actionClass with these args but don't start it */
    public StsAction newAction(Class actionClass, Object[] args)
    {
        try
        {
            if (actionClass == null) return null;
            this.actionClass = actionClass;

            Object[] arguments = getArgs(args);
            Class[] argTypes = getArgTypes(arguments);

            Constructor constructor = actionClass.getDeclaredConstructor(argTypes);
            if (constructor == null)
            {
                StsException.systemError("constructor not found for action " + actionClass.getName());
                return null;
            }
            StsAction newAction = (StsAction) constructor.newInstance(arguments);
            newAction.setActionButtonText();
            return newAction;
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage(actionClass.getName() + " action failed. " + e.toString());
            StsException.outputException("StsActionManager.newAction() failed.", e, StsException.WARNING);
            return null;
        }
    }

    /** create instance of actionClass and start it */
    public StsAction startAction(Class actionClass)
    {
        StsAction newAction = newAction(actionClass);
        return startAction(newAction, null);
    }

    /** create instance of actionClass and start it */
    public void startAction(Class actionClass, StsAction lastAction)
    {
        StsAction newAction = newAction(actionClass);
        startAction(newAction, lastAction);
    }

    /** create instance of actionClass with these args and start it */
    public StsAction startAction(Class actionClass, Object[] args)
    {
        try
        {
            StsAction action = newAction(actionClass, args);
            return checkStartAction(action);
        }
        catch (Exception e)
        {
            StsMessageFiles.logMessage(actionClass.getName() + " action failed. " + e.toString());
            StsException.outputException("StsActionManager.newAction() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public StsAction checkStartAction(StsAction newAction)
    {
        if (newAction != null && newAction.checkStartAction())
            return startAction(newAction, getCurrentAction());
        else
            return null;
    }

    /**
     * lastAction may have been completed; if this is a repeat of same action, classInitialize it as a repeatAction.
     * If the new action is passive (canInterrupt==true), simply push it onto the activeList stack.
     * Otherwise if the currentAction can be terminated ok or the user wants to terminate it,
     * then end it and start the newAction.
     */
    public StsAction startAction(StsAction newAction, StsAction lastAction)
    {
        try
        {
            if (newAction == null) return null;
            if (actionDebug) StsException.systemDebug(this, "startAction ", " currentAction: " + newAction.getName());
            //            if(lastAction != null && (!newAction.canInterrupt() || !lastAction.canBeInterrupted())) return null;

            logUsage(START, newAction);

            if (lastAction != null && lastAction.getClass() == newAction.getClass())
                newAction.initializeRepeatAction(lastAction);
            else
            {
                StsAction currentAction = getCurrentAction();

                if (currentAction != null)
                {
                    if (currentAction.isRepeatable())
                        currentAction.setIsRepeatable(false);
                    if (!newAction.canInterrupt()) endCurrentAction();
                }
                addNewAction(newAction);
            }

            // check if this action wants to open a transaction (generally for the life of the action)
            // or transactions will be issued for each object created by the action
            // most actions do the former.  If you wish the latter, set openTransaction to false
            // in the constructor of the action
            //           if(newAction.hasOpenTransaction())
            model.initializeTransaction(newAction.getClass().getName());

            newAction.checkAddToolbar();

            //            checkStartAuxiliaryWindowAction(newAction);

            if (!newAction.checkPrerequisites())
                return null;

            if (newAction instanceof Runnable)
            {
                Thread thread = new Thread((Runnable) newAction);
                thread.start();
            }
            else if (!newAction.start())
            {
                abortCurrentAction();
                removeCurrentAction();
            }

            if (newAction != null)
                newAction.setStarted(true);

            return newAction;
        }
        catch (Exception e)
        {
            StsMessageFiles.errorMessage("StsActionManager.startAction() failed for action: " + StsToolkit.getSimpleClassname(newAction));
            StsException.outputException("StsActionManager.startAction() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public void addNewAction(StsAction action)
    {
        if (actionDebug) StsException.systemDebug(this, "addNewAction ", " Current action: " + action.getName());
        removeExtraActionHack(action);
        actionList.add(action);
        action.setActionButtonText();
    }

    /** logs usage for currentAction. State is is a string: Start, End, Abort. */
    private void logUsage(int state, StsAction currentAction)
    {
        String actionWorkflowStepName = getWorkflowStepName(currentAction);
        switch (state)
        {
            case START:
                // Log with old module and message
                // Set new module and message
                if (debug)
                    System.out.println("START log - Module: " + Main.usageModule + " Message: " + Main.usageMessage);
                Main.logUsage(currentAction.getClass().getName(), actionWorkflowStepName);
                if (debug)
                    System.out.println("START set to - Module: " + Main.usageModule + " Message: " + Main.usageMessage);
                break;
            case END:
                // Log with old module and message
                if (debug) System.out.println("END - Module: " + Main.usageModule + " Message: " + Main.usageMessage);
                Main.logUsage();
                Main.logReset();
                if (debug)
                    System.out.println("END set to - Module: " + Main.usageModule + " Message: " + Main.usageMessage);
                break;
            case ABORT:
                // Log with old module
                if (debug) System.out.println("ABORT - Module: " + Main.usageModule + " Message: " + Main.usageMessage);
                Main.logUsage();
                Main.logReset();
                if (debug)
                    System.out.println("ABORT set to - Module: " + Main.usageModule + " Message: " + Main.usageMessage);
                break;
        }
    }

    /** logs usage at fixed intervals for this action */
    public void logUsageTimer(StsAction currentAction)
    {
        String actionWorkflowStepName = getWorkflowStepName(currentAction);
        //        System.out.println("Left Mouse - Module: " + Main.usageModule + " Message: " + Main.usageMessage);
        Main.logUsageTimer(currentAction.getClass().getName(), actionWorkflowStepName);
        //        System.out.println("Left Mouse set to - Module: " + Main.usageModule + " Message: " + Main.usageMessage);
    }

    public void endAction(Class actionClass)
    {
        if (actionClass.isInstance(getCurrentAction()))
            terminateCurrentAction();
    }

    public void endAction(StsAction action)
    {
        if (getCurrentAction() == action) endCurrentAction();
    }

    public void endAction()
    {
        terminateCurrentAction();
    }

    private boolean endCurrentActionOK()
    {
        StsAction currentAction = getCurrentAction();
        if (currentAction == null) return true;

        // if the currentAction can be interrupted,
        // leave it on the actionList stack; it will be reactivated when newAction completes
        if (currentAction.canBeInterrupted())
        {
            //            endCurrentAction();
            return true;
        }

        // if this action is not running an open transaction, don't bother the user
        //        if(!currentAction.hasOpenTransaction()) return true;

        boolean ok = StsYesNoDialog.questionValue(currentAction.getModel().win3d, StsToolkit.getSimpleClassname(currentAction) +
            " has not completed, Do you want to end it?");
        if (!ok) return false;

        endCurrentAction();
        return true;
    }

    public void removeCurrentAction()
    {
        StsAction currentAction = getCurrentAction();
        if (currentAction == null) return;
        currentAction.setActionButtonTextToEnd();
        removeExtraActionHack(currentAction);
    }

    public void removeAction(StsAction action)
    {
        if(action == null) return;
        if (actionDebug) StsException.systemDebug(this, "removeAction ", "Current action: " + action.getName());
        removeExtraActionHack(action);
        //        actionList.remove(action);
        action.setActionButtonTextToEnd();
    }

    /**
     * sombody is putting the same action on the action list twice;
     * until we figure it out, we need to remove all actions of this type when
     * we start a newAction as well as when we end an action.
     *
     * @param action
     */
    private void removeExtraActionHack(StsAction action)
    {
        Class actionClass = action.getClass();
        int nCurrentActions = actionList.size();
        int nRemoved = 0;
        for (int n = nCurrentActions - 1; n >= 0; n--)
        {
            StsAction listAction = (StsAction) actionList.get(n);
            if (listAction.getClass() == actionClass)
            {
                actionList.remove(listAction);
                nRemoved++;
            }
        }
        if (nRemoved > 1)
            StsException.systemDebug(this, "removeExtraActionHack ", "Current action: " + action.getName());
    }

    public boolean abortCurrentTransaction()
    {
        // abort the current transaction
        try
        {
            if (model.currentTransaction == null) return true;
            return model.currentTransaction.abort();
        }
        catch (Exception e)
        {
            StsException.outputException("StsActionManager.abortCurrentTransaction() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    public void undoAction()
    {
        abortCurrentTransaction();
        repaint();
    }

    public void repaint()
    {
        model.win3dDisplayAll();
    }

    public void abortCurrentAction()
    {
        StsAction currentAction = getCurrentAction();

        if (currentAction == null) return;

        try
        {
            model.currentTransaction.abort();
            logUsage(ABORT, currentAction);
            removeCurrentAction();

        }
        catch (Exception e)
        {
            StsException.systemError("StsActionManager.abortCurrentAction() failed.");
        }
    }

}
