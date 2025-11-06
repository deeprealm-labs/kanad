/**
 * Professional Card Component - IBM-Inspired Design
 *
 * Features:
 * - Sharp rectangular design (no rounded corners)
 * - IBM-inspired color palette
 * - Optional header, body, footer sections
 * - Minimal shadow (IBM style)
 * - Hover effect support
 */

import React from 'react';

export interface CardProps {
  children: React.ReactNode;
  className?: string;
  hover?: boolean;
}

export interface CardHeaderProps {
  children: React.ReactNode;
  className?: string;
}

export interface CardBodyProps {
  children: React.ReactNode;
  className?: string;
}

export interface CardFooterProps {
  children: React.ReactNode;
  className?: string;
}

export interface CardTitleProps {
  children: React.ReactNode;
  className?: string;
}

// Main Card Component
const Card: React.FC<CardProps> & {
  Header: React.FC<CardHeaderProps>;
  Body: React.FC<CardBodyProps>;
  Footer: React.FC<CardFooterProps>;
  Title: React.FC<CardTitleProps>;
} = ({ children, className = '', hover = false }) => {
  const cardClasses = `
    card
    ${hover ? 'cursor-pointer' : ''}
    ${className}
  `.trim().replace(/\s+/g, ' ');

  return <div className={cardClasses}>{children}</div>;
};

// Card Header
const CardHeader: React.FC<CardHeaderProps> = ({ children, className = '' }) => {
  return <div className={`card-header ${className}`}>{children}</div>;
};

// Card Body
const CardBody: React.FC<CardBodyProps> = ({ children, className = '' }) => {
  return <div className={`card-body ${className}`}>{children}</div>;
};

// Card Footer
const CardFooter: React.FC<CardFooterProps> = ({ children, className = '' }) => {
  return <div className={`card-footer ${className}`}>{children}</div>;
};

// Card Title
const CardTitle: React.FC<CardTitleProps> = ({ children, className = '' }) => {
  return <h3 className={`card-title ${className}`}>{children}</h3>;
};

// Attach subcomponents
Card.Header = CardHeader;
Card.Body = CardBody;
Card.Footer = CardFooter;
Card.Title = CardTitle;

export default Card;
