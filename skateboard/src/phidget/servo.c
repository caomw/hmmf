#include <stdio.h>
#include <phidget21.h>

int AttachHandler(CPhidgetHandle SERV, void *userptr)
{
	int serialNo;
	const char *name;

	CPhidget_getDeviceName (SERV, &name);
	CPhidget_getSerialNumber(SERV, &serialNo);
	printf("%s %10d attached!\n", name, serialNo);

	return 0;
}

int DetachHandler(CPhidgetHandle SERV, void *userptr)
{
	int serialNo;
	const char *name;

	CPhidget_getDeviceName (SERV, &name);
	CPhidget_getSerialNumber(SERV, &serialNo);
	printf("%s %10d detached!\n", name, serialNo);

	return 0;
}

int ErrorHandler(CPhidgetHandle SERV, void *userptr, int ErrorCode, const char *Description)
{
	printf("Error handled. %d - %s\n", ErrorCode, Description);
	return 0;
}

int PositionChangeHandler(CPhidgetServoHandle SERV, void *usrptr, int Index, double Value)
{
	printf("Motor: %d > Current Position: %f\n", Index, Value);
	return 0;
}

void attach_handlers(CPhidgetHandle SERV)
{
	//Set the handlers to be run when the device is plugged in or opened from software, unplugged or closed from software, or generates an error.
	CPhidget_set_OnAttach_Handler( SERV, AttachHandler, NULL);
	CPhidget_set_OnDetach_Handler( SERV, DetachHandler, NULL);
	CPhidget_set_OnError_Handler( SERV, ErrorHandler, NULL);
	
    //Registers a callback that will run when the motor position is changed.
	//Requires the handle for the Phidget, the function that will be called, and an arbitrary pointer that will be supplied to the callback function (may be NULL).
	CPhidgetServo_set_OnPositionChange_Handler( (CPhidgetServoHandle) SERV, PositionChangeHandler, NULL);
}

int open_servo( CPhidgetHandle SERV)
{
	int result;
	const char *err;
	//open the servo for device connections
	CPhidget_open( SERV, -1);

	//get the program to wait for an servo device to be attached
	printf("Waiting for Servo controller to be attached....");
	if((result = CPhidget_waitForAttachment( SERV, 10000)))
	{
		CPhidget_getErrorDescription(result, &err);
		printf("Problem waiting for attachment: %s\n", err);
		return 0;
	}

}

//Display the properties of the attached phidget to the screen.  We will be displaying the name, serial number and version of the attached device.
int display_properties(CPhidgetServoHandle phid)
{
	int serialNo, version, numMotors;
	const char* ptr;

	CPhidget_getDeviceType((CPhidgetHandle)phid, &ptr);
	CPhidget_getSerialNumber((CPhidgetHandle)phid, &serialNo);
	CPhidget_getDeviceVersion((CPhidgetHandle)phid, &version);

	CPhidgetServo_getMotorCount (phid, &numMotors);

	printf("%s\n", ptr);
	printf("Serial Number: %10d\nVersion: %8d\n# Motors: %d\n", serialNo, version, numMotors);

	return 0;
}

