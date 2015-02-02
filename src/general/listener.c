/*

Copyright (c) 2003 - 2005 by:
Leif Laaksonen , Centre for Scientific Computing , ESPOO, FINLAND
Confidential unpublished property of Leif Laaksonen
All rights reserved

Coded by:
Eero Häkkinen
*/

#include "maindefs.h"

#include <assert.h>
#include <stdlib.h>

#include "gomlistener.h"
#include "printmsg.h"

#include "stdafx.h"

/***************************************************************************/
gom_Listener* gomp_AddListener(
    gom_ListenerList *listPointer,
    gom_ListenerFunc  callback,
    void             *callbackData )
/***************************************************************************/
{
    gom_Listener *listener;
    listener = malloc(sizeof(*listener));
    if ( ! listener ) {
        gomp_PrintERROR("Unable to allocate memory for a listener.");
        return(NULL);
    }
    assert(callback != NULL);
    /* Store the data. */
    listener->prev         = NULL;
    listener->next         = *listPointer;
    if ( listener->next )
        listener->next->prev = listener;
    listener->callback     = callback;
    listener->callbackData = callbackData;
    *listPointer           = listener;
    return(listener);
}
/***************************************************************************/
static void CancelListener(
    gom_ListenerList *listPointer,
    gom_Listener     *listener )
/***************************************************************************/
{
    if ( listener->prev )
        listener->prev->next = listener->next;
    else {
        assert(*listPointer == listener);
        *listPointer = listener->next;
    }
    if ( listener->next )
        listener->next->prev = listener->prev;
}
/***************************************************************************/
int gomp_CancelListener(
    gom_ListenerList *listPointer,
    gom_Listener     *listener )
/***************************************************************************/
{
    if ( listener ) {
        CancelListener(listPointer, listener);
        free(listener);
        return(1);
    }
    return(0);
}
/***************************************************************************/
int gomp_CancelListenersByFunc(
    gom_ListenerList *listPointer,
    gom_ListenerFunc  callback )
/***************************************************************************/
{
    gom_Listener *listener = *listPointer;
    gom_Listener *next;
    int count = 0;

    assert(callback != NULL);

    while ( listener ) {
        next = listener->next;
        if ( listener->callback == callback ) {
            /* Remove the listener from the list. */
            CancelListener(listPointer, listener);
            free(listener);
            count++;
        }
        listener = next;
    }
    
    return count;
}
/***************************************************************************/
int gomp_CallListeners(
    gom_ListenerList      *listPointer,
    gom_ListenerHelperFunc helper,
    void                  *helperData )
/***************************************************************************/
{
    gom_Listener *listener = *listPointer;
    gom_Listener  local;
    int failed = 0, result;

    if ( ! listener )
        /* We have no listeners to call. */
        return(0);

    /* Prevent gomp_CancelListenersByFunc to remove &local. */
    local.callback = NULL;

    while ( listener ) {
        /* Add &local to the listener list after listener.
        ** No one except us can cancel &local so after the call of
        ** CallListener local.next is still valid.
        **/
        local.prev       = listener;
        local.next       = listener->next;
        local.prev->next = &local;
        if ( local.next )
            local.next->prev = &local;

        result = helper(
            helperData,
            listener->callback,
            listener->callbackData );
        if ( result != 0 ) {
            if ( result > 0 )
                /* An error occured. */
                failed = 1;
            else if ( local.prev == listener ) {
                /* Remove the listener from the list. */
                CancelListener(listPointer, listener);
                free(listener);
            }
        }

        /* Get the next listener. */
        listener = local.next;

        /* Remove the temporary listener from the list. */
        CancelListener(listPointer, &local);
    }

    return failed;
}
/***************************************************************************/
static int CallSimpleListener(
    void            *helperData,
    gom_ListenerFunc callback,
    void            *callbackData )
/***************************************************************************/
{
    return ((gom_SimpleListenerFunc)callback)(callbackData);
}
/***************************************************************************/
int gomp_CallSimpleListeners(
    gom_ListenerList *listPointer )
/***************************************************************************/
{
    if ( ! *listPointer )
        /* We have no listeners to call. */
        return(0);
    return gomp_CallListeners(listPointer,CallSimpleListener,NULL);
}
/***************************************************************************/
/*
** Implement update data listeners.
**/
static int gomp_CallUpdateDataListeners( void );
GOM_IMPLEMENT_SIMPLE_GOM_LISTENER(UpdateData)
/***************************************************************************/
int gomp_UpdateData(void)
/***************************************************************************/
{
    return gomp_CallUpdateDataListeners();
}
/***************************************************************************/
