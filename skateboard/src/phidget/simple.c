#include "servo.c"


int servo_simple()
{
	double curr_pos = 0.0;

	// servo handle
	CPhidgetServoHandle servo = 0;

	//create object
	CPhidgetServo_create(&servo);

    // attach all handlers
    attach_handlers( (CPhidgetHandle) servo);
    
    // open servo
    open_servo( (CPhidgetHandle) servo);

	//Display the properties of the attached servo device
	display_properties(servo);

	//read servo event data
	printf("Reading.....\n");

	// Motor at index 0
	//display current motor position
	CPhidgetServo_getPosition(servo, 0, &curr_pos);
	printf("Motor: 0 > Current Position: %f\n", curr_pos);

	//keep displaying servo event data until user input is read
	printf("Press any key to continue\n");
	getchar();

	//change the motor position
	//valid range is -22 to 232
	//we'll set it to a few random positions to move it around

	//Step 1: Position 10.00
	printf("Move to position 10.00. Press any key to Continue\n");
	getchar();
	CPhidgetServo_setPosition (servo, 0, 10.00);

	printf("Move to position 100.00. Press any key to Continue\n");
	getchar();
	CPhidgetServo_setPosition (servo, 0, 100.00);

	//Step 5: Position 200.00
	printf("Move to position 200.00. Press any key to Continue\n");
	getchar();
	CPhidgetServo_setPosition (servo, 0, 200.00);

	//Step 6: Position 20.00
	printf("Move to position 20.00. Press any key to Continue\n");
	getchar();

	CPhidgetServo_setPosition (servo, 0, 20.00);

	//Step 7: Diseangage
	printf("Disengage. Press any key to Continue\n");
	getchar();
	CPhidgetServo_setEngaged (servo, 0, 0);

	printf("Press any key to end\n");
	getchar();

	//since user input has been read, this is a signal to terminate the program so we will close the phidget and delete the object we created
	printf("Closing...\n");
	CPhidget_close((CPhidgetHandle)servo);
	CPhidget_delete((CPhidgetHandle)servo);

	//all done, exit
	return 0;
}

int main(int argc, char* argv[])
{
	servo_simple();
	return 0;
}

