// SPDX-FileCopyrightText: 2026 Starsworn (fork of Erin Catto's Box3D)
// SPDX-License-Identifier: MIT

// A gear joint: two bodies' spins about their own z axes held in a ratio, as meshed teeth hold them. One velocity
// constraint on the angular velocities only:
//   radiusA * dot(wA, axisA) + radiusB * dot(wB, axisB) = carrierSpeed
// where carrierSpeed is (radiusA + radiusB) times the spin of what both gears are mounted on (the game sets it each
// step: gears on a turning chassis). Its torque is limited (teeth that strip past it). No recording ops: a replay
// skips gear joints.

#include "body.h"
#include "core.h"
#include "joint.h"
#include "math_internal.h"
#include "physics_world.h"
#include "solver.h"
#include "solver_set.h"

// needed for dll export
#include "box3d/box3d.h"

void b3GearJoint_SetCarrierSpeed( b3JointId jointId, float speed )
{
	B3_ASSERT( b3IsValidFloat( speed ) );
	b3JointSim* base = b3GetJointSimCheckType( jointId, b3_gearJoint );
	base->gearJoint.carrierSpeed = speed;
}

float b3GearJoint_GetCarrierSpeed( b3JointId jointId )
{
	b3JointSim* base = b3GetJointSimCheckType( jointId, b3_gearJoint );
	return base->gearJoint.carrierSpeed;
}

void b3GearJoint_SetMaxTorque( b3JointId jointId, float torque )
{
	B3_ASSERT( b3IsValidFloat( torque ) && torque >= 0.0f );
	b3JointSim* base = b3GetJointSimCheckType( jointId, b3_gearJoint );
	base->gearJoint.maxTorque = torque;
}

float b3GearJoint_GetMeshTorque( b3JointId jointId )
{
	b3World* world = b3GetWorld( jointId.world0 );
	b3JointSim* base = b3GetJointSimCheckType( jointId, b3_gearJoint );
	return world->inv_h * base->gearJoint.impulse * base->gearJoint.radiusA;
}

b3Vec3 b3GetGearJointTorque( b3World* world, b3JointSim* base )
{
	b3GearJoint* joint = &base->gearJoint;
	return b3MulSV( world->inv_h * joint->impulse * joint->radiusA, joint->axisA );
}

void b3PrepareGearJoint( b3JointSim* base, b3StepContext* context )
{
	B3_ASSERT( base->type == b3_gearJoint );

	b3World* world = context->world;

	b3Body* bodyA = b3Array_Get( world->bodies, base->bodyIdA );
	b3Body* bodyB = b3Array_Get( world->bodies, base->bodyIdB );

	B3_ASSERT( bodyA->setIndex == b3_awakeSet || bodyB->setIndex == b3_awakeSet );
	b3SolverSet* setA = b3Array_Get( world->solverSets, bodyA->setIndex );
	b3SolverSet* setB = b3Array_Get( world->solverSets, bodyB->setIndex );

	int localIndexA = bodyA->localIndex;
	int localIndexB = bodyB->localIndex;

	b3BodySim* bodySimA = b3Array_Get( setA->bodySims, localIndexA );
	b3BodySim* bodySimB = b3Array_Get( setB->bodySims, localIndexB );

	base->invMassA = bodySimA->invMass;
	base->invMassB = bodySimB->invMass;
	base->invIA = bodySimA->invInertiaWorld;
	base->invIB = bodySimB->invInertiaWorld;

	b3GearJoint* joint = &base->gearJoint;
	joint->indexA = bodyA->setIndex == b3_awakeSet ? localIndexA : B3_NULL_INDEX;
	joint->indexB = bodyB->setIndex == b3_awakeSet ? localIndexB : B3_NULL_INDEX;

	b3Quat quatA = b3MulQuat( bodySimA->transform.q, base->localFrameA.q );
	b3Quat quatB = b3MulQuat( bodySimB->transform.q, base->localFrameB.q );
	joint->axisA = b3RotateVector( quatA, b3Vec3_axisZ );
	joint->axisB = b3RotateVector( quatB, b3Vec3_axisZ );

	b3Vec3 jA = b3MulSV( joint->radiusA, joint->axisA );
	b3Vec3 jB = b3MulSV( joint->radiusB, joint->axisB );
	float k = b3Dot( jA, b3MulMV( base->invIA, jA ) ) + b3Dot( jB, b3MulMV( base->invIB, jB ) );
	joint->mass = k > 0.0f ? 1.0f / k : 0.0f;

	if ( context->enableWarmStarting == false )
	{
		joint->impulse = 0.0f;
	}
}

static void b3ApplyGearImpulse( b3JointSim* base, b3StepContext* context, float impulse )
{
	b3GearJoint* joint = &base->gearJoint;
	b3BodyState dummyState = b3_identityBodyState;
	b3BodyState* stateA = joint->indexA == B3_NULL_INDEX ? &dummyState : context->states + joint->indexA;
	b3BodyState* stateB = joint->indexB == B3_NULL_INDEX ? &dummyState : context->states + joint->indexB;

	b3Vec3 jA = b3MulSV( joint->radiusA * impulse, joint->axisA );
	b3Vec3 jB = b3MulSV( joint->radiusB * impulse, joint->axisB );

	if ( stateA->flags & b3_dynamicFlag )
	{
		stateA->angularVelocity = b3Add( stateA->angularVelocity, b3MulMV( base->invIA, jA ) );
	}
	if ( stateB->flags & b3_dynamicFlag )
	{
		stateB->angularVelocity = b3Add( stateB->angularVelocity, b3MulMV( base->invIB, jB ) );
	}
}

void b3WarmStartGearJoint( b3JointSim* base, b3StepContext* context )
{
	B3_ASSERT( base->type == b3_gearJoint );
	b3ApplyGearImpulse( base, context, base->gearJoint.impulse );
}

void b3SolveGearJoint( b3JointSim* base, b3StepContext* context )
{
	B3_ASSERT( base->type == b3_gearJoint );
	b3GearJoint* joint = &base->gearJoint;

	b3BodyState dummyState = b3_identityBodyState;
	b3BodyState* stateA = joint->indexA == B3_NULL_INDEX ? &dummyState : context->states + joint->indexA;
	b3BodyState* stateB = joint->indexB == B3_NULL_INDEX ? &dummyState : context->states + joint->indexB;

	float cdot = joint->radiusA * b3Dot( stateA->angularVelocity, joint->axisA ) +
				 joint->radiusB * b3Dot( stateB->angularVelocity, joint->axisB ) - joint->carrierSpeed;

	float maxImpulse = context->h * joint->maxTorque / b3MaxFloat( joint->radiusA, 1.0e-3f );
	float oldImpulse = joint->impulse;
	joint->impulse = b3ClampFloat( oldImpulse - joint->mass * cdot, -maxImpulse, maxImpulse );
	b3ApplyGearImpulse( base, context, joint->impulse - oldImpulse );
}

void b3DrawGearJoint( b3DebugDraw* draw, b3JointSim* base, b3WorldTransform transformA, b3WorldTransform transformB, float scale )
{
	float length = 0.1f * scale;
	b3WorldTransform frameA = b3MulWorldTransforms( transformA, base->localFrameA );
	b3WorldTransform frameB = b3MulWorldTransforms( transformB, base->localFrameB );
	draw->DrawSegmentFcn( frameA.p, frameB.p, b3_colorGold, draw->context );
	draw->DrawSegmentFcn( frameA.p, b3OffsetPos( frameA.p, b3MulSV( length, b3RotateVector( frameA.q, b3Vec3_axisZ ) ) ), b3_colorGreen, draw->context );
	draw->DrawSegmentFcn( frameB.p, b3OffsetPos( frameB.p, b3MulSV( length, b3RotateVector( frameB.q, b3Vec3_axisZ ) ) ), b3_colorBlue, draw->context );
}
